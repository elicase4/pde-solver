namespace pdesolver::application::heateq::problem {

	template<typename Backend, typename HeatEqBundle>
	HeatProblem<Backend, HeatEqBundle>::HeatProblem(const application::heateq::config::HeatConfig& config, mesh::Mesh mesh) : config_(config), mesh_(std::move(mesh)), topoDOF_(mesh_, config_.discretization.dofOrdering), sourceForms_(config_.source.expression) {

		// conductivity model
		if (config_.conductivity.type != config::ConductivityConfig::Type::Constant) {
			throw std::runtime_error("HeatProblem: only constant conductivity is supported so far");
		}
		conductivityModel_.conductivity = config_.conductivity.value;

		// boundary conditions
		Index numFluxBCs = config_.boundaryConditions.size();
		fluxForms_.reserve(numFluxBCs);
		
		for (const auto& bcCfg : config_.boundaryConditions) {
			for (auto form : bcCfg.forms) {
				
				if (form == application::heateq::config::BoundaryConditionConfig::Form::ValueBC) {
					
					using DirichletT = typename HeatEqBundle::template DirichletBC<SourceCallableT>;

					auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<DirichletT>>(
						new fem::boundary::BoundaryCondition<DirichletT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Essential}, DirichletT{bcCfg.expression}});
					essentialBCs_.registerBC<DirichletT>(bc);
				
				} else if (form == config::BoundaryConditionConfig::Form::FluxBC) {

					using FluxFunctionT = typename HeatEqBundle::template FluxBC<FluxCallableT>;

					fluxForms_.push_back(std::make_unique<FluxFormsT>(bcCfg.expression));

					auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<FluxFunctionT>>(
						new fem::boundary::BoundaryCondition<FluxFunctionT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Natural}, FluxFunctionT{bcCfg.expression}});

					naturalBCs_.template registerBC<FluxFunctionT>(bc, *fluxForms_.back(), defaultModelBdy_);
				
				}
			
			}
		}

		// build topoDOF constrains and linear system containers
		topoDOF_.template buildConstraints<typename HeatEqBundle::Basis>(essentialBCs_);
		
		if (config_.solver.linear->operatorType == solver::config::LinearSolverConfig::OperatorType::CSR) {
			K_ = std::make_unique<MatrixT>(fem::assembly::Assembler<Backend>::template createMatrix<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));
		}

		F_ = std::make_unique<VectorT>(fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));

		U_ = std::make_unique<VectorT>(fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));
		U_->zero();

		// TODO: branching and advanced handling for linear, non linear solvers, steady/transient
		if (config_.solver.linear->operatorType == solver::config::LinearSolverConfig::OperatorType::CSR) {
		    CSROperatorT op(*K_);
		    linearSolverRunner_ = solver::linear::makeLinearSolverRunner<CSROperatorT, VectorT>(op, topoDOF_.numFreeDOFs(), *config_.solver.linear, "heateq");
		} else if (config_.solver.linear->operatorType == solver::config::LinearSolverConfig::OperatorType::FEM) {
		    FEMOperatorT op(assembler_, mesh_, topoDOF_, Real(0), conductivityModel_, matrixForms_);
		    linearSolverRunner_ = solver::linear::makeLinearSolverRunner<FEMOperatorT, VectorT>(op, topoDOF_.numFreeDOFs(), *config_.solver.linear, "heateq");
		} else {
			throw std::runtime_error("HeatProblem: unsupported operator");
		}
	
	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::assembleSystem(Real time) {

		if (config_.solver.linear->operatorType == solver::config::LinearSolverConfig::OperatorType::CSR) {
		    fem::assembly::Assembler<Backend>::template assembleMatrix<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::ConstantConductivityModel, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, conductivityModel_, matrixForms_, *U_, *K_);
		}

		F_->zero();
		fem::assembly::Assembler<Backend>::template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::DefaultModel, SourceFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, defaultModel_, sourceForms_, *U_, *F_);

		bcApplicator_.template applyNaturalBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPBdy, typename HeatEqBundle::QuadratureBoundaryType>(mesh_, topoDOF_, naturalBCs_, time, *F_);

		bcApplicator_.template applyEssentialBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::ConstantConductivityModel, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, essentialBCs_, time, conductivityModel_, matrixForms_, *F_);
	
	}

	template<typename Backend, typename HeatEqBundle>
	bool HeatProblem<Backend, HeatEqBundle>::solveLinear() {

		linalg::solver::SolverReport<VectorT> report;
		bool converged = linearSolverRunner_->solve(*F_, *U_, report);
		return converged;
	
	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::writeOutput(Index step) const {

		if (!config_.output.vtk || (step % config_.output.writeFrequency != 0)) return;
		
		const std::string filename = config_.output.directory + "/" + config_.output.prefix + "_" + std::to_string(step) + ".vtk";
		io::FieldIO::writeVTK<HeatEqBundle::NumDOFs>(mesh_, topoDOF_, essentialBCs_, Real(0), U_->data(), {"T"}, filename);
		
	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::writeLog() const {

	}

} // namespace pdesolver::application::heateq::problem
