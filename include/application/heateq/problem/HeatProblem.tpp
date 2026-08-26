namespace pdesolver::application::heateq::problem {

	template<typename Backend, typename HeatEqBundle>
	HeatProblem<Backend, HeatEqBundle>::HeatProblem(const application::heateq::config::HeatConfig& config, mesh::Mesh mesh, typename HeatEqBundle::Basis basis, typename HeatEqBundle::QuadratureVolumeType quadratureVolume, typename HeatEqBundle::QuadratureBoundaryType quadratureBoundary) : config_(config), solverInstance_(solver::resolveSolverInstance(config_.solver)), mesh_(std::move(mesh)), topoDOF_(mesh_, config_.discretization.dofOrdering), driverLogger_(solver::logging::makeDriverLogger(config_.logging.driver, "heateq")), evalEleTemplate_(std::move(basis)), quadratureVolume_(std::move(quadratureVolume)), quadratureBoundary_(std::move(quadratureBoundary)), sourceForms_(config_.source.expression) {

		// solver instance
		if (solver::isTransient(solverInstance_.mode)) {
			throw std::runtime_error("HeatProblem: transient driver support is not yet implemented");
		}

		// conductivity model
		if (config_.conductivity.type == config::ConductivityConfig::Type::Constant) {
			conductivityModel_ = typename HeatEqBundle::ConstantConductivityModel{config_.conductivity.value};
		} else if (config_.conductivity.type == config::ConductivityConfig::Type::Anisotropic) {
			conductivityModel_ = typename HeatEqBundle::AnisotropicConductivityModel{config_.conductivity.tensor};
		} else {
			throw std::runtime_error("HeatProblem: unsupported conductivity type");
		}

		// boundary conditions
		Index numFluxBCs = config_.boundaryConditions.size();
		fluxForms_.reserve(numFluxBCs);
		
		for (const auto& bcCfg : config_.boundaryConditions) {
			for (auto form : bcCfg.forms) {

				if (form == application::heateq::config::BoundaryConditionConfig::Form::ValueBC) {

					if (bcCfg.mode == config::BoundaryConditionConfig::Mode::File) {
						throw std::runtime_error("HeatProblem: file-based boundary conditions are not yet implemented");
					}

					using DirichletT = typename HeatEqBundle::template DirichletBC<SourceCallableT>;

					auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<DirichletT>>(new fem::boundary::BoundaryCondition<DirichletT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Essential}, DirichletT{bcCfg.expression}});

					essentialBCs_.registerBC<DirichletT>(bc);

				} else if (form == config::BoundaryConditionConfig::Form::FluxBC) {

					if (bcCfg.mode == config::BoundaryConditionConfig::Mode::File) {
						throw std::runtime_error("HeatProblem: file-based boundary conditions are not yet implemented");
					}

					using FluxFunctionT = typename HeatEqBundle::template FluxBC<FluxCallableT>;

					if (bcCfg.fluxExpression.size() != HeatEqBundle::SpatialDim) {
						throw std::runtime_error("HeatProblem: flux boundary 'expression' must have SpatialDim components");
					}

					fluxForms_.push_back(std::make_unique<FluxFormsT>(bcCfg.fluxExpression));

					auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<FluxFunctionT>>(new fem::boundary::BoundaryCondition<FluxFunctionT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Natural}, FluxFunctionT{bcCfg.fluxExpression}});

					naturalBCs_.template registerBC<FluxFunctionT>(bc, *fluxForms_.back(), defaultModelBdy_);

				}

			}
		}

		// build topoDOF constrains and linear system containers
		topoDOF_.buildConstraints(evalEleTemplate_.basis(), essentialBCs_);
		
		if (solverInstance_.linear->operatorType == solver::config::LinearSolverConfig::OperatorType::CSR) {
			K_ = std::make_unique<MatrixT>(fem::assembly::Assembler<Backend>::template createMatrix<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));
		}

		F_ = std::make_unique<VectorT>(fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));

		U_ = std::make_unique<VectorT>(fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));
		U_->zero();

		// parse intital condition
		if (config_.initialCondition.type == config::InitialConditionConfig::Type::Expression) {

			utils::expression::ScalarExpression icExpr(config_.initialCondition.expression);

			for (Index nodeID = 0; nodeID < mesh_.data.numNodes; ++nodeID) {

				Real coords[3] = {Real(0), Real(0), Real(0)};
				const Real* nodeCoordPtr = mesh_.getNodeCoord(nodeID);
				for (Index d = 0; d < HeatEqBundle::SpatialDim; ++d) coords[d] = nodeCoordPtr[d];

				Real icVal[HeatEqBundle::NumDOFs];
				icExpr(Real(0), coords, icVal);

				for (Index c = 0; c < HeatEqBundle::NumDOFs; ++c) {

					Index tdof = topoDOF_.getNodeDOF(nodeID, c);
					if (topoDOF_.isConstrained(tdof)) continue;

					Index adof = topoDOF_.toAlgebraic(tdof);
					U_->data()[adof] = icVal[c];

				}

			}

		} else {
			throw std::runtime_error("HeatProblem: file-based initial conditions are not yet implemented");
		}

		if (solverInstance_.linear->operatorType == solver::config::LinearSolverConfig::OperatorType::CSR) {
		    
			CSROperatorT op(*K_);
			linearSolverRunner_ = solver::linear::makeLinearSolverRunner<CSROperatorT, VectorT>(op, topoDOF_.numFreeDOFs(), *solverInstance_.linear, config_.logging.solver, "Heat Equation", std::vector<std::string>{"T"});
		
		} else if (solverInstance_.linear->operatorType == solver::config::LinearSolverConfig::OperatorType::FEM) {

			std::visit([&](auto& model) {

				using ConductivityModelT = std::decay_t<decltype(model)>;
				using FEMOperatorT = FEMOperatorFor<ConductivityModelT>;

				FEMOperatorT op(assembler_, mesh_, topoDOF_, Real(0), model, matrixForms_, evalEleTemplate_, quadratureVolume_);
				linearSolverRunner_ = solver::linear::makeLinearSolverRunner<FEMOperatorT, VectorT>(op, topoDOF_.numFreeDOFs(), *solverInstance_.linear, config_.logging.solver, "Heat Equation", std::vector<std::string>{"T"});

			}, conductivityModel_);

		} else {
			
			throw std::runtime_error("HeatProblem: unsupported operator");
		
		}
	
	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::assembleSystem(Real time) {

		if (solverInstance_.linear->operatorType == solver::config::LinearSolverConfig::OperatorType::CSR) {
			std::visit([&](auto& model) {
				using ConductivityModelT = std::decay_t<decltype(model)>;
				fem::assembly::Assembler<Backend>::template assembleMatrix<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, ConductivityModelT, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, model, matrixForms_, evalEleTemplate_, quadratureVolume_, *U_, *K_);
			}, conductivityModel_);
		}

		F_->zero();
		fem::assembly::Assembler<Backend>::template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::DefaultModel, SourceFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, defaultModel_, sourceForms_, evalEleTemplate_, quadratureVolume_, *U_, *F_);

		bcApplicator_.template applyNaturalBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPBdy, typename HeatEqBundle::QuadratureBoundaryType>(mesh_, topoDOF_, naturalBCs_, time, evalEleTemplate_, quadratureBoundary_, *F_);

		std::visit([&](auto& model) {
			using ConductivityModelT = std::decay_t<decltype(model)>;
			bcApplicator_.template applyEssentialBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, ConductivityModelT, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, essentialBCs_, time, model, matrixForms_, evalEleTemplate_, quadratureVolume_, *F_);
		}, conductivityModel_);
	
	}

	template<typename Backend, typename HeatEqBundle>
	bool HeatProblem<Backend, HeatEqBundle>::solveLinear() {

		linalg::solver::SolverReport<VectorT> report;
		bool converged = linearSolverRunner_->solve(*F_, *U_, report);
		return converged;
	
	}

	template<typename Backend, typename HeatEqBundle>
	Real HeatProblem<Backend, HeatEqBundle>::residualNorm() const {

		// PSEUDOCODE -- unreachable until Newton/Picard exist (see resolveSolverInstance() note
		// in the constructor). K_/F_/U_ reflect whatever was last assembled by assembleSystem();
		// this must NOT re-assemble itself, callers are expected to call assemble() first.
		//
		// r = F_ - K_ * U_        (raw allocation TBD -- reuse a scratch VectorT member rather
		//                          than allocating here, this runs every nonlinear iteration)
		// return ||r||_2

		return Real(0);

	}

	template<typename Backend, typename HeatEqBundle>
	bool HeatProblem<Backend, HeatEqBundle>::solveLinearStep() {

		// PSEUDOCODE -- unreachable until Newton/Picard exist.
		//
		// Newton (solverInstance_.nonlinear->type == Newton):
		//   assemble the EXACT tangent J(U_) -- NOT the same as K_ once a NonlinearTangentForm-
		//   conforming form exists for a genuinely U-dependent equation (K_ alone is missing the
		//   dK/dU * U term). Heat conduction has no such form yet (constant conductivity), so
		//   this path has nothing correct to fall back to -- that's exactly why Newton stays
		//   unreachable rather than silently reusing K_ as if it were the true tangent.
		//
		// Picard (solverInstance_.nonlinear->type == Picard):
		//   J(U_) == K_ (already assembled, U-dependent coefficients frozen at current U_) --
		//   this path COULD reuse the existing linearSolverRunner_ machinery directly:
		//     r = F_ - K_ * U_
		//     solve K_ * deltaU = r   via linearSolverRunner_
		//     U_ += deltaU
		//     return whether that linear solve converged
		//
		// Either way: this is a correction solve (deltaU), NOT solveLinear()'s direct solve for
		// U_ -- do not collapse the two.

		return false;

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::writeOutput(Index step) const {

		if (!config_.output.vtk || (step % config_.output.writeFrequency != 0)) return;

		const std::string filename = config_.output.directory + "/" + config_.output.prefix + "_" + std::to_string(step) + ".vtk";
		io::FieldIO::writeVTK<HeatEqBundle::NumDOFs>(mesh_, topoDOF_, essentialBCs_, Real(0), U_->data(), {"T"}, filename);

		const char* ordStr = (topoDOF_.ordering() == fem::dof::DOFOrdering::Interleaved) ? "Interleaved" : "Block";
		driverLogger_.event("wrote '" + filename + "' - " + std::to_string(HeatEqBundle::NumDOFs) + " field(s), " + std::to_string(mesh_.data.numNodes) + " nodes, " + ordStr + " ordering");

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::writeLog() const {

	}

} // namespace pdesolver::application::heateq::problem
