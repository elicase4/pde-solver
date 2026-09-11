namespace pdesolver::application::heateq::problem {

	template<typename Backend, typename HeatEqBundle>
	HeatProblem<Backend, HeatEqBundle>::HeatProblem(const application::heateq::config::HeatConfig& config, mesh::Mesh mesh, typename HeatEqBundle::Basis basis, typename HeatEqBundle::QuadratureVolumeType quadratureVolume, typename HeatEqBundle::QuadratureBoundaryType quadratureBoundary) : config_(config), solverInstance_(solver::resolveSolverInstance(config_.solver)), mesh_(std::move(mesh)), topoDOF_(mesh_, config_.discretization.dofOrdering), driverLogger_(solver::logging::makeDriverLogger(config_.logging.driver, "heateq")), evalEleTemplate_(std::move(basis)), quadratureVolume_(std::move(quadratureVolume)), quadratureBoundary_(std::move(quadratureBoundary)) {

		// conductivity model
		if (config_.conductivity.type == config::ConductivityConfig::Type::Constant) {
			conductivityModel_.setConstant(config_.conductivity.value);
			conductivityModelBdy_.setConstant(config_.conductivity.value);
		} else if (config_.conductivity.type == config::ConductivityConfig::Type::Anisotropic) {
			conductivityModel_.setAnisotropic(config_.conductivity.tensor);
			conductivityModelBdy_.setAnisotropic(config_.conductivity.tensor);
		} else {
			throw std::runtime_error("HeatProblem: unsupported conductivity type");
		}

		const bool transient = solver::isTransient(solverInstance_.mode);

		// transient setup
		if (transient) {

			if (!config_.density.has_value() || !config_.specificHeat.has_value()) {
				throw std::runtime_error("HeatProblem: 'materials.density' and 'materials.specific_heat' are required for a transient driver");
			}

			if (!solverInstance_.timestepper) {
				throw std::runtime_error("HeatProblem: solver.timestepper config is required for a transient driver");
			}

			densityModel_.value = config_.density->value;
			specificHeatModel_.value = config_.specificHeat->value;

			materialModel_ = typename HeatEqBundle::MaterialModel(conductivityModel_, densityModel_, specificHeatModel_);
			massMaterialModel_ = typename HeatEqBundle::MassMaterialModel(densityModel_, specificHeatModel_);

			setDt(solverInstance_.timestepper->stepSize.dt);

		}

		// monitors
		for (const auto& monitorCfg : config_.monitors) {

			const std::string unit = fem::quantity::unitFor<typename HeatEqBundle::HeatFluxIntegrand>(monitorCfg.reduction, HeatEqBundle::SpatialDim - 1);
			const std::vector<std::string> columns{"tick", "time", monitorCfg.name + " [" + unit + "]"};

			if (monitorCfg.reduction == fem::quantity::Reduction::Integral) {

				MonitorCombinationIntegralT combination;
				for (const auto& term : monitorCfg.terms) {
					monitorRegistryIntegral_.registerTag(term.boundary);
					combination.addTerm(term.boundary, term.coefficient);
				}
				monitorOutputsIntegral_.push_back(MonitorOutput<MonitorCombinationIntegralT>{monitorCfg.name, std::move(combination), monitorCfg.output.console, unit, utils::logging::CsvWriter(monitorCfg.output.file, columns)});

			} else {

				MonitorCombinationAverageT combination;
				for (const auto& term : monitorCfg.terms) {
					monitorRegistryAverage_.registerTag(term.boundary);
					combination.addTerm(term.boundary, term.coefficient);
				}
				monitorOutputsAverage_.push_back(MonitorOutput<MonitorCombinationAverageT>{monitorCfg.name, std::move(combination), monitorCfg.output.console, unit, utils::logging::CsvWriter(monitorCfg.output.file, columns)});

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

		// constraints & linear system containers
		topoDOF_.buildConstraints(evalEleTemplate_.basis(), essentialBCs_);

		if (solverInstance_.linear->operatorType == OpType::CSR) {
			K_ = std::make_unique<MatrixT>(fem::assembly::Assembler<Backend>::template createMatrix<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));
		}

		F_ = std::make_unique<VectorT>(fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));
		U_ = std::make_unique<VectorT>(fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));
		U_->zero();

		if (transient) {
			U_prev_ = std::make_unique<VectorT>(fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));
			massTimesUprev_ = std::make_unique<VectorT>(fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_));
		}

		// initial condition
		loadNodalField(config_.initialCondition.read, *U_);
		if (transient) {
			linalg::operations::copy(*U_, *U_prev_);
		}

		// configure the linear operator
		const OpType opType = solverInstance_.linear->operatorType;

		if (opType == OpType::CSR) {

			if (transient){
				assembleTransientOperatorMatrix();
			}
			makeLinearRunner(CSROperatorT(*K_));

		} else if (opType == OpType::FEM) {

			if (transient) {
				makeLinearRunner(TransientMatrixFreeOperatorT(assembler_, mesh_, topoDOF_, Real(0), materialModel_, transientOperatorForms_, evalEleTemplate_, quadratureVolume_));
			} else {
				makeLinearRunner(SteadyMatrixFreeOperatorT(assembler_, mesh_, topoDOF_, Real(0), conductivityModel_, stiffnessForms_, evalEleTemplate_, quadratureVolume_));
			}

		} else {
			throw std::runtime_error("HeatProblem: unsupported operator");
		}

	}

	// assembly helpers
	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::loadNodalField(const solver::config::NodalFieldReadConfig& cfg, VectorT& target) const {

		auto scatter = [&](Index nodeID, const Real* val) {
			for (Index c = 0; c < HeatEqBundle::NumDOFs; ++c) {
				Index tdof = topoDOF_.getNodeDOF(nodeID, c);
				if (topoDOF_.isConstrained(tdof)) continue;
				target.data()[topoDOF_.toAlgebraic(tdof)] = val[c];
			}
		};

		if (cfg.mode == solver::config::NodalFieldReadConfig::Mode::Expression) {

			utils::expression::ScalarExpression expr(cfg.expression);

			for (Index nodeID = 0; nodeID < mesh_.data.numNodes; ++nodeID) {
				Real coords[3] = {Real(0), Real(0), Real(0)};
				const Real* p = mesh_.getNodeCoord(nodeID);
				for (Index d = 0; d < HeatEqBundle::SpatialDim; ++d) coords[d] = p[d];
				Real val[HeatEqBundle::NumDOFs];
				expr(Real(0), coords, val);
				scatter(nodeID, val);
			}

		} else {

			NodalScalarSourceT src(mesh_, cfg.file);

			for (Index nodeID = 0; nodeID < mesh_.data.numNodes; ++nodeID) {
				Real val[HeatEqBundle::NumDOFs];
				src.eval(nodeID, val);
				scatter(nodeID, val);
			}

		}

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::assembleStiffnessMatrix(Real time) {

		fem::assembly::Assembler<Backend>::template assembleMatrix<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::ConductivityModel, StiffnessFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, conductivityModel_, stiffnessForms_, evalEleTemplate_, quadratureVolume_, *U_, *K_);

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::assembleTransientOperatorMatrix() {

		fem::assembly::Assembler<Backend>::template assembleMatrix<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::MaterialModel, TransientOperatorFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, Real(0), materialModel_, transientOperatorForms_, evalEleTemplate_, quadratureVolume_, *U_, *K_);

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::assembleLoad(Real time) {

		F_->zero();

		if (config_.source.read.mode == solver::config::NodalFieldReadConfig::Mode::Expression) {
			fem::assembly::Assembler<Backend>::template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::DefaultModel, ExpressionSourceFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, defaultModel_, *sourceForms_, evalEleTemplate_, quadratureVolume_, *U_, *F_);
		} else {
			fem::assembly::Assembler<Backend>::template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::DefaultModel, NodalSourceFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, time, defaultModel_, *nodalSourceForms_, evalEleTemplate_, quadratureVolume_, *U_, *F_);
		}

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::addTransientMassTerm(Real time) {

		massTimesUprev_->zero();
		fem::assembly::Assembler<Backend>::template assembleVector<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::MassMaterialModel, MassFormsT, typename HeatEqBundle::QuadratureVolumeType, fem::assembly::GatherMode::Full>(mesh_, topoDOF_, time - dt_, massMaterialModel_, massForms_, evalEleTemplate_, quadratureVolume_, *U_prev_, *massTimesUprev_, &essentialBCs_);
		linalg::operations::axpy(transientMassCoefficient(), *massTimesUprev_, *F_);

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::applyNatural(Real time) {

		bcApplicator_.template applyNaturalBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPBdy, typename HeatEqBundle::QuadratureBoundaryType>(mesh_, topoDOF_, naturalBCs_, time, evalEleTemplate_, quadratureBoundary_, *F_);

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::applyEssential(Real time) {

		if (solver::isTransient(solverInstance_.mode)) {
			bcApplicator_.template applyEssentialBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::MaterialModel, TransientOperatorFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, essentialBCs_, time, materialModel_, transientOperatorForms_, evalEleTemplate_, quadratureVolume_, *F_);
		} else {
			bcApplicator_.template applyEssentialBCs<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::ConductivityModel, StiffnessFormsT, typename HeatEqBundle::QuadratureVolumeType>(mesh_, topoDOF_, essentialBCs_, time, conductivityModel_, stiffnessForms_, evalEleTemplate_, quadratureVolume_, *F_);
		}

	}

	template<typename Backend, typename HeatEqBundle>
	template<typename OperatorT>
	void HeatProblem<Backend, HeatEqBundle>::makeLinearRunner(const OperatorT& op) {

		linearSolverRunner_ = solver::linear::makeLinearSolverRunner<OperatorT, VectorT>(op, topoDOF_.numFreeDOFs(), *solverInstance_.linear, config_.logging.solver, "Heat Equation", std::vector<std::string>{"T"});

	}

	// public surface
	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::assembleSystem(Real time) {

		const bool transient = solver::isTransient(solverInstance_.mode);

		if (!transient && solverInstance_.linear->operatorType == OpType::CSR) {
			assembleStiffnessMatrix(time);
		}

		assembleLoad(time);
		if (transient){
			addTransientMassTerm(time);
		}

		applyNatural(time);
		applyEssential(time);

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::advanceTimestep() {

		// step-to-step state roll -- driven by the timestepper via HeatStage::advance()
		linalg::operations::copy(*U_, *U_prev_);

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::setDt(Real dt) {

		if (dt == dt_) return; // common case: the active step-size policy kept dt unchanged

		dt_ = dt;
		transientOperatorForms_ = TransientOperatorFormsT(typename HeatEqBundle::DiffusionForm{}, typename HeatEqBundle::MassFormOverDt{transientMassCoefficient()});

		if (K_) {
			assembleTransientOperatorMatrix();
		}

	}

	template<typename Backend, typename HeatEqBundle>
	bool HeatProblem<Backend, HeatEqBundle>::solveLinear() {

		linalg::solver::SolverReport<VectorT> report;
		return linearSolverRunner_->solve(*F_, *U_, report);

	}

	template<typename Backend, typename HeatEqBundle>
	Real HeatProblem<Backend, HeatEqBundle>::residualNorm() const {

		// PSEUDOCODE -- unreachable until Newton/Picard exist. K_/F_/U_ reflect the last
		// assembleSystem(); this must NOT re-assemble. Callers assemble() first.
		//   r = F_ - K_ * U_   (reuse a scratch VectorT member, not a fresh allocation)
		//   return ||r||_2

		return Real(0);

	}

	template<typename Backend, typename HeatEqBundle>
	bool HeatProblem<Backend, HeatEqBundle>::solveLinearStep() {

		// PSEUDOCODE -- unreachable until Newton/Picard exist.
		//
		// Newton: assemble the EXACT tangent J(U_) -- not K_, which is missing the dK/dU * U term.
		// Heat conduction has no U-dependent form yet, so there's nothing correct to fall back to.
		//
		// Picard: J(U_) == K_ (coefficients frozen at U_) -- could reuse linearSolverRunner_:
		//   r = F_ - K_ * U_ ; solve K_ * deltaU = r ; U_ += deltaU ; return converged
		//
		// Either way this is a correction solve (deltaU), NOT solveLinear()'s direct solve.

		return false;

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::writeOutput(Index step) const {

		if (!config_.output.vtk || (step % config_.output.writeFrequency != 0)) return;

		const std::string filename = config_.output.directory + "/" + config_.output.prefix + "_" + std::to_string(step) + ".vtk";
		io::fieldio::FieldIO::writeVTK<HeatEqBundle::NumDOFs>(mesh_, topoDOF_, essentialBCs_, Real(0), U_->data(), {"T [K]"}, filename);

		const char* ordStr = (topoDOF_.ordering() == fem::dof::DOFOrdering::Interleaved) ? "Interleaved" : "Block";
		driverLogger_.event("wrote '" + filename + "' - " + std::to_string(HeatEqBundle::NumDOFs) + " field(s), " + std::to_string(mesh_.data.numNodes) + " nodes, " + ordStr + " ordering");

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::writeLog() const {

	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::evaluateMonitors(Index tick, Real time) const {

		if (monitorOutputsIntegral_.empty() && monitorOutputsAverage_.empty()) return;

		if (!monitorOutputsIntegral_.empty()) {
			fem::quantity::QuantityEvaluator<Backend>::template evaluateBoundaryRegistry<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPBdy, typename HeatEqBundle::ConductivityModelBdy, MonitorQuantitiesIntegralT, typename HeatEqBundle::QuadratureBoundaryType>(mesh_, topoDOF_, essentialBCs_, time, conductivityModelBdy_, monitorQuantitiesIntegral_, evalEleTemplate_, quadratureBoundary_, *U_, monitorRegistryIntegral_);
		}

		if (!monitorOutputsAverage_.empty()) {
			fem::quantity::QuantityEvaluator<Backend>::template evaluateBoundaryRegistry<HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPBdy, typename HeatEqBundle::ConductivityModelBdy, MonitorQuantitiesAverageT, typename HeatEqBundle::QuadratureBoundaryType>(mesh_, topoDOF_, essentialBCs_, time, conductivityModelBdy_, monitorQuantitiesAverage_, evalEleTemplate_, quadratureBoundary_, *U_, monitorRegistryAverage_);
		}

		Real value[HeatEqBundle::HeatFluxIntegrand::NumComponents];

		for (const auto& out : monitorOutputsIntegral_) {
			out.combination.evaluate(monitorRegistryIntegral_, value);
			if (out.toConsole) driverLogger_.event("monitor '" + out.name + "' = " + std::to_string(value[0]) + " " + out.unit);
			out.csv.writeRow({static_cast<Real>(tick), time, value[0]});
		}

		for (const auto& out : monitorOutputsAverage_) {
			out.combination.evaluate(monitorRegistryAverage_, value);
			if (out.toConsole) driverLogger_.event("monitor '" + out.name + "' = " + std::to_string(value[0]) + " " + out.unit);
			out.csv.writeRow({static_cast<Real>(tick), time, value[0]});
		}

	}

} // namespace pdesolver::application::heateq::problem
