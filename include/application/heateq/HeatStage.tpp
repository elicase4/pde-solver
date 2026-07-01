namespace pdesolver::application::heateq {

	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::HeatStage(const HeatConfig& config) : config_(config) {}

	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	void HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::initialize() {

		// Load mesh from binary PMSH file.
		if (config_.mesh.file.empty()) {
			throw std::runtime_error("HeatStage::initialize: mesh.file is empty — no generator path implemented yet");
		}
		
		io::MeshIO::readBinary(mesh_, config_.mesh.file);

		if (!mesh_.isValid()) {
			throw std::runtime_error("HeatStage::initialize: loaded mesh failed validity check");
		}

		fem::dof::DOFOrdering dofOrder = config_.discretization.blockDOFOrdering ? fem::dof::DOFOrdering::Block : fem::dof::DOFOrdering::Interleaved;

		// Build topological DOF manager
		topoDOF_ = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh_, dofOrder);

		// build source function
		auto srcExpr = std::make_shared<utils::expression::ScalarExpression>(condif_.source.expression);
		auto sourceFunc_ = [srcExpr](Real t, const Real* x, Real* out) {
			(*srcExpr)(t, x, out);
		};

		// Register BCs
		for (const auto& bcCfg : config_.boundaryConditions) {
			
			if (bcCfg.type == BoundaryConditionConfig::Type::Value) {
				auto expr = std::make_shared<utils::expression::ScalarExpression>(bcCfg.expression);
				auto bcFunction = [expr](Real t, const Real* x, Real* out) { (*expr)(t,x,out); };
				using BCType = typename HeatEqBundle::template DirichletBC<decltype(bcFunction)>;
				using BCCat = fem::boundary::BCCategory::Essential;
				BCType bc(bcFunction);
				fem::boundary::BoundaryCondition<BcType> wrappedBC(bcCfg.boundaryID, BCCat, bc);
				bcRegistry_.registerBC(wrappedBC);
			} else if (bcCfg.type == BoundaryConditionConfig::Type::Flux) {
				auto expr = std::make_shared<utils::expression::VectorExpression>(bcCfg.expression);
				auto bcFunction = [expr](Real t, const Real* x, Real* out) { (*expr)(t,x,out); };
				using BCType = typename HeatEqBundle::template FluxBC<decltype(bcFunction)>;
				using BCCat = fem::boundary::BCCategory::Natural;
				BCType bc(bcFunction);
				fem::boundary::BoundaryCondition<BcType> wrappedBC(bcCfg.boundaryID, BCCat, bc);
				bcRegistry_.registerBC(wrappedBC);
			}


		}

		topoDOF_->buildConstraints<BasisType>(bcRegistry_);

		// Allocate system
		if (!config_.linearSolver.matrixFree) {
			K_ = assembler_.template createMatrix<HeatEqBundle::NumDOFs>(mesh_, *topoDOF_);
		}
		
		F_ = assembler_.template createVector<HeatEqBundle::NumDOFs>(mesh_, *topoDOF_);
		U_ = assembler_.template createVector<HeatEqBundle::NumDOFs>(mesh_, *topoDOF_);

	}

	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	template<>
	void HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::assemble<true>() {

		const Real t = 0.0; // steady

		// Assemble stiffness matrix
		typename HeatEqBundle::DiffusionForm diffusionForm;
		fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms(diffusionForm);

		assembler_.template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::ConstantConductivityModel, decltype(operatorForms), QuadratureVolumeType>(mesh_, *topoDOF_, t, conductivityModel_, operatorForms, U_, O_);
		
		// Assemble RHS vector
		HeatEqBundle::SourceFunction<decltype(sourceFunc)> sourceFunction(sourceFunc);
		HeatEqBundle::SourceForm<decltype(sourceFunc)> sourceForm(sourceFunction);
		fem::form::FormRegistry<HeatEqBundle::SourceForm<decltype(f)>> rhsForms(sorurceForm);
		
		assembler_.template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::DefaultModel, decltype(rhsForms), QuadratureVolumeType>(mesh_, *topoDOF_, t, defaultModel_, rhsForms, U_, F_);

		// Apply essential BCs
		bcApplicator_.template applyEssentialBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, decltype(operatorForms), typename HeatEqBundle::ConstantConductivityModel, QuadratureVolumeType>(mesh_, *topoDOF_, bcRegistry_, t, conductivityModel_, operatorForms, F_);
	
		// Apply natural BCs
		bcApplicator_.template applyNaturalBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, decltype(fluxForms_), QuadratureBoundaryType>(mesh_, *topoDOF_, bcRegistry_, t, fluxForms_, F_);
	
	}

	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	template<>
	void HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::assemble<false>() {

		const Real t = 0.0; // steady

		// Assemble stiffness matrix
		typename HeatEqBundle::DiffusionForm diffusionForm;

		assembler_.template assembleMatrix<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::ConstantConductivityModel, typename HeatEqBundle::DiffusionForm, QuadratureVolumeType>(mesh_, *topoDOF_, t, conductivityModel_, diffusionForm, U_, K_);
		
		// Assemble RHS vector
		HeatEqBundle::SourceFunction<decltype(sourceFunc)> sourceFunction(sourceFunc);
		HeatEqBundle::SourceForm<decltype(sourceFunc)> sourceForm(sourceFunction);
		
		assembler_.template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::DefaultModel, typename HeatEqBundle::SourceForm<decltype()>, QuadratureVolumeType>(mesh_, *topoDOF_, t, defaultModel_, sourceForm, U_, F_);

		// Apply essential BCs
		bcApplicator_.template applyEssentialBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::DiffusionForm, typename HeatEqBundle::ConstantConductivityModel, QuadratureVolumeType>(mesh_, *topoDOF_, bcRegistry_, t, conductivityModel_, diffusionForm, F_);
	
		// Apply natural BCs
		bcApplicator_.template applyNaturalBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::FluxForm<>, QuadratureBoundaryType>(mesh_, *topoDOF_, bcRegistry_, t, fluxForm, F_);
	
	}
	
	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	template<>
	bool HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::solve<true>() {

		using Workspace = linalg::solver::iterative::cg::Workspace<VectorT>;
		using Report = linalg::solver::SolverReport<VectorT>;
		using Preconditioner = linalg::solver::preconditioner::Identity<VectorT>;
		using Logger = utils::logging::ConsoleLogger;
		using Config = linalg::solver::iterative::cg::Config<VectorT>;
		using Solver = linalg::solver::iterative::cg::Solver<Operator, VectorT, Preconditioner, Logger>;

		Operator op(K_);
		Workspace W(topoDOF_->numFreeDOFs());
		Report report;
		Preconditioner M;
		Logger logger("CG", "theta");

		CGConfig cgCfg{config_.linearSolver.tolerance, linalg::solver::iterative::cg::ToleranceType::Relative, config_.linearSolver.maxIterations};

		CGSolver solver(cgCfg);
		solver.solve(report, logger, W, M, op, F_, U_);

		return report.converged;
	}
	
	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	template<>
	bool HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::solve<true>() {

		using Workspace = linalg::solver::iterative::cg::Workspace<VectorT>;
		using Report = linalg::solver::SolverReport<VectorT>;
		using Preconditioner = linalg::solver::preconditioner::Identity<VectorT>;
		using Logger = utils::logging::ConsoleLogger;
		using Config = linalg::solver::iterative::cg::Config<VectorT>;
		using Solver = linalg::solver::iterative::cg::Solver<Operator, VectorT, Preconditioner, Logger>;

		Operator op(O_);
		Workspace W(topoDOF_->numFreeDOFs());
		Report report;
		Preconditioner M;
		Logger logger("CG", "theta");

		CGConfig cgCfg{config_.linearSolver.tolerance, linalg::solver::iterative::cg::ToleranceType::Relative, config_.linearSolver.maxIterations};

		CGSolver solver(cgCfg);
		solver.solve(report, logger, W, M, op, F_, U_);

		return report.converged;
	}

	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	void HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::finalize() {

		if (!config_.output.vtk) {
			return;
		}

		const Real t = 0.0;

		std::filesystem::create_directories(config_.output.directory);

		const auto path = std::filesystem::path(config_.output.directory) / (config_.output.prefix + ".vtk");

		io::FieldIO::writeVTK(mesh_, *topoDOF_, bcRegistry_, t, U_.data(), {"T"}, path.string());
	}

} // namespace pdesolver::application::heateq
