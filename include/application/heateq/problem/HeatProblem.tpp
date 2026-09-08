namespace pdesolver::application::heateq::problem {

	template<typename Backend, typename HeatEqBundle>
	HeatProblem<Backend, HeatEqBundle>::HeatProblem(const application::heateq::config::HeatConfig& config, mesh::Mesh mesh, typename HeatEqBundle::Basis basis, typename HeatEqBundle::QuadratureVolumeType quadratureVolume, typename HeatEqBundle::QuadratureBoundaryType quadratureBoundary) : config_(config), solverInstance_(solver::resolveSolverInstance(config_.solver)), mesh_(std::move(mesh)), topoDOF_(mesh_, config_.discretization.dofOrdering), driverLogger_(solver::logging::makeDriverLogger(config_.logging.driver, "heateq")), evalEleTemplate_(std::move(basis)), quadratureVolume_(std::move(quadratureVolume)), quadratureBoundary_(std::move(quadratureBoundary)) {

		// solver instance
		if (solver::isTransient(solverInstance_.mode)) {
			throw std::runtime_error("HeatProblem: transient driver support is not yet implemented");
		}

		// conductivity model
		if (config_.conductivity.type == config::ConductivityConfig::Type::Constant) {
			conductivityModel_ = typename HeatEqBundle::ConstantConductivityModel{config_.conductivity.value};
			conductivityModelBdy_ = typename HeatEqBundle::ConstantConductivityModelBdy{config_.conductivity.value};
		} else if (config_.conductivity.type == config::ConductivityConfig::Type::Anisotropic) {
			conductivityModel_ = typename HeatEqBundle::AnisotropicConductivityModel{config_.conductivity.tensor};
			conductivityModelBdy_ = typename HeatEqBundle::AnisotropicConductivityModelBdy{config_.conductivity.tensor};
		} else {
			throw std::runtime_error("HeatProblem: unsupported conductivity type");
		}

		// monitors
		for (const auto& monitorCfg : config_.monitors) {

			if (monitorCfg.reduction == fem::quantity::Reduction::Integral) {

				MonitorCombinationIntegralT combination;
				for (const auto& term : monitorCfg.terms) {
					monitorRegistryIntegral_.registerTag(term.boundary);
					combination.addTerm(term.boundary, term.coefficient);
				}
				monitorCombinationsIntegral_.emplace_back(monitorCfg.name, std::move(combination));

			} else {

				MonitorCombinationAverageT combination;
				for (const auto& term : monitorCfg.terms) {
					monitorRegistryAverage_.registerTag(term.boundary);
					combination.addTerm(term.boundary, term.coefficient);
				}
				monitorCombinationsAverage_.emplace_back(monitorCfg.name, std::move(combination));

			}

		}

		// source
		if (config_.source.read.mode == solver::config::NodalFieldReadConfig::Mode::Expression) {
			sourceForms_.emplace(config_.source.read.expression);
		} else {
			nodalSourceForms_.emplace(NodalScalarSourceT(mesh_, config_.source.read.file));
		}

		// boundary conditions
		Index numFluxBCs = config_.boundaryConditions.size();
		expressionFluxForms_.reserve(numFluxBCs);
		nodalFluxForms_.reserve(numFluxBCs);
		
		for (const auto& bcCfg : config_.boundaryConditions) {
			
			for (auto form : bcCfg.forms) {

				if (form == application::heateq::config::BoundaryConditionConfig::Form::ValueBC) {

					if (bcCfg.mode == solver::config::NodalFieldReadConfig::Mode::File) {

						auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<DirichletNodalT>>(new fem::boundary::BoundaryCondition<DirichletNodalT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Essential}, DirichletNodalT{mesh_, bcCfg.file}});

						essentialBCs_.registerBC<DirichletNodalT>(bc);

					} else {

						auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<DirichletExpressionT>>(new fem::boundary::BoundaryCondition<DirichletExpressionT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Essential}, DirichletExpressionT{bcCfg.expression}});

						essentialBCs_.registerBC<DirichletExpressionT>(bc);

					}

				} else if (form == config::BoundaryConditionConfig::Form::FluxBC) {

					if (bcCfg.mode == solver::config::NodalFieldReadConfig::Mode::File) {

						NodalFluxSourceT nodalFluxSource(mesh_, bcCfg.file);
						nodalFluxForms_.push_back(std::make_unique<NodalFluxFormsT>(nodalFluxSource));

						auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<FluxFunctionNodalT>>(new fem::boundary::BoundaryCondition<FluxFunctionNodalT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Natural}, FluxFunctionNodalT{mesh_, bcCfg.file}});

						naturalBCs_.template registerBC<FluxFunctionNodalT>(bc, *nodalFluxForms_.back(), defaultModelBdy_);

					} else {

						if (bcCfg.fluxExpression.size() != HeatEqBundle::SpatialDim) {
							throw std::runtime_error("HeatProblem: flux boundary 'expression' must have SpatialDim components");
						}

						expressionFluxForms_.push_back(std::make_unique<ExpressionFluxFormsT>(bcCfg.fluxExpression));

						auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<FluxFunctionExpressionT>>(new fem::boundary::BoundaryCondition<FluxFunctionExpressionT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Natural}, FluxFunctionExpressionT{bcCfg.fluxExpression}});

						naturalBCs_.template registerBC<FluxFunctionExpressionT>(bc, *expressionFluxForms_.back(), defaultModelBdy_);

					}

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

		// initial condition
		if (config_.initialCondition.read.mode == solver::config::NodalFieldReadConfig::Mode::Expression) {

			utils::expression::ScalarExpression icExpr(config_.initialCondition.read.expression);

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

			NodalScalarSourceT icSource(mesh_, config_.initialCondition.read.file);

			for (Index nodeID = 0; nodeID < mesh_.data.numNodes; ++nodeID) {

				Real icVal[HeatEqBundle::NumDOFs];
				icSource.eval(nodeID, icVal);

				for (Index c = 0; c < HeatEqBundle::NumDOFs; ++c) {

					Index tdof = topoDOF_.getNodeDOF(nodeID, c);
					if (topoDOF_.isConstrained(tdof)) continue;

					Index adof = topoDOF_.toAlgebraic(tdof);
					U_->data()[adof] = icVal[c];

				}

			}

		}

		// configure linear operator
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

		// assemble K if needed
		if (solverInstance_.linear->operatorType == solver::config::LinearSolverConfig::OperatorType::CSR) {
			std::visit([&](auto& model) {
				using ConductivityModelT = std::decay_t<decltype(model)>;
				fem::assembly::Assembler<Backend>::template assembleMatrix<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, ConductivityModelT, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, model, matrixForms_, evalEleTemplate_, quadratureVolume_, *U_, *K_);
			}, conductivityModel_);
		}

		// assemble F
		F_->zero();
		if (config_.source.read.mode == solver::config::NodalFieldReadConfig::Mode::Expression) {
			fem::assembly::Assembler<Backend>::template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::DefaultModel, ExpressionSourceFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, defaultModel_, *sourceForms_, evalEleTemplate_, quadratureVolume_, *U_, *F_);
		} else {
			fem::assembly::Assembler<Backend>::template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::DefaultModel, NodalSourceFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, defaultModel_, *nodalSourceForms_, evalEleTemplate_, quadratureVolume_, *U_, *F_);
		}

		// apply natural BCs
		bcApplicator_.template applyNaturalBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPBdy, typename HeatEqBundle::QuadratureBoundaryType>(mesh_, topoDOF_, naturalBCs_, time, evalEleTemplate_, quadratureBoundary_, *F_);

		// apply essential BCs
		std::visit([&](auto& model) {
			using ConductivityModelT = std::decay_t<decltype(model)>;
			bcApplicator_.template applyEssentialBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, ConductivityModelT, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, essentialBCs_, time, model, matrixForms_, evalEleTemplate_, quadratureVolume_, *F_);
		}, conductivityModel_);
	
	}

	template<typename Backend, typename HeatEqBundle>
	bool HeatProblem<Backend, HeatEqBundle>::solveLinear() {

		// run linear solver
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
		io::fieldio::FieldIO::writeVTK<HeatEqBundle::NumDOFs>(mesh_, topoDOF_, essentialBCs_, Real(0), U_->data(), {"T"}, filename);

		const char* ordStr = (topoDOF_.ordering() == fem::dof::DOFOrdering::Interleaved) ? "Interleaved" : "Block";
		driverLogger_.event("wrote '" + filename + "' - " + std::to_string(HeatEqBundle::NumDOFs) + " field(s), " + std::to_string(mesh_.data.numNodes) + " nodes, " + ordStr + " ordering");

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::writeLog() const {

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::evaluateMonitors() const {

		if (monitorCombinationsIntegral_.empty() && monitorCombinationsAverage_.empty()) return;

		std::visit([&](auto& model) {

			using ConductivityModelT = std::decay_t<decltype(model)>;

			if (!monitorCombinationsIntegral_.empty()) {
				fem::quantity::QuantityEvaluator<Backend>::template evaluateBoundaryRegistry<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPBdy, ConductivityModelT, MonitorQuantitiesIntegralT, typename HeatEqBundle::QuadratureBoundaryType>(mesh_, topoDOF_, essentialBCs_, Real(0), model, monitorQuantitiesIntegral_, evalEleTemplate_, quadratureBoundary_, *U_, monitorRegistryIntegral_);
			}

			if (!monitorCombinationsAverage_.empty()) {
				fem::quantity::QuantityEvaluator<Backend>::template evaluateBoundaryRegistry<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPBdy, ConductivityModelT, MonitorQuantitiesAverageT, typename HeatEqBundle::QuadratureBoundaryType>(mesh_, topoDOF_, essentialBCs_, Real(0), model, monitorQuantitiesAverage_, evalEleTemplate_, quadratureBoundary_, *U_, monitorRegistryAverage_);
			}

		}, conductivityModelBdy_);

		Real value[HeatEqBundle::HeatFluxIntegrand::NumComponents];

		for (const auto& [name, combination] : monitorCombinationsIntegral_) {
			combination.evaluate(monitorRegistryIntegral_, value);
			driverLogger_.event("monitor '" + name + "' = " + std::to_string(value[0]));
		}

		for (const auto& [name, combination] : monitorCombinationsAverage_) {
			combination.evaluate(monitorRegistryAverage_, value);
			driverLogger_.event("monitor '" + name + "' = " + std::to_string(value[0]));
		}

	}

} // namespace pdesolver::application::heateq::problem
