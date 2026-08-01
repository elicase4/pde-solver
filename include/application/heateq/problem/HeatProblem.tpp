namespace pdesolver::application::heateq::problem {

	template<typename Backend, typename HeatEqBundle>
	HeatProblem<Backend, HeatEqBundle>::HeatProblem(const config::HeatConfig& config, mesh::Mesh mesh) : config_(config), mesh_(std::move(mesh)), topoDOF_(mesh_, config_.discretization.dofOrdering) {

		// assembler_, bcApplicator_, conductivityModel_, defaultModel_, defaultModelBdy_,
		// matrixForms_ all default-construct fine (all stateless except
		// conductivityModel_, which needs its .conductivity set in the body below).
		
		// --- 1. conductivity model -------------------------------------------------
		// V1: only ConductivityConfig::Type::Constant is wired up (matches
		// constantConductivityModel.conductivity = 1.0 in the integration tests).
		// Throw for Anisotropic/Default until a real model backs them -- see the
		// open extensibility question about a swappable conductivity-model axis on
		// HeatEqBundle discussed with the user; not resolved yet, deliberately
		// hardcoded to Constant for now.
		//
		//   if (config_.conductivity.type != config::ConductivityConfig::Type::Constant)
		//       throw std::runtime_error("HeatProblem: only constant conductivity is supported so far");
		//   conductivityModel_.conductivity = config_.conductivity.value;

		// --- 2. source term ----------------------------------------------------------
		// Mirrors the integration tests' `sourceFunction`/`sourceForm`/`rhsForms` triple,
		// just built from config instead of a lambda:
		//
		//   using SourceFunctionT = typename HeatEqBundle::template SourceFunction<SourceCallableT>;
		//   using SourceFormT = typename HeatEqBundle::template SourceForm<SourceCallableT>;
		//   SourceFunctionT sourceFunction{SourceCallableT(config_.source.expression)};
		//   sourceForms_ = SourceFormsT{SourceFormT{sourceFunction}};

		// --- 3. boundary conditions --------------------------------------------------
		// fluxForms_ MUST be reserve()'d before any push_back/emplace_back --
		// NaturalBoundaryOperator stores a reference into whichever element it was
		// registered with, so a mid-loop reallocation would leave every
		// already-registered natural BC dangling.
		//
		//   Index numFluxBCs = count of entries in config_.boundaryConditions whose
		//                      forms include config::BoundaryConditionConfig::Form::FluxBC;
		//   fluxForms_.reserve(numFluxBCs);
		//
		//   for (const auto& bcCfg : config_.boundaryConditions) {
		//       for (auto form : bcCfg.forms) {
		//           if (form == config::BoundaryConditionConfig::Form::ValueBC) {
		//               // essential (Dirichlet) BC -- mirrors bc0/bc1/bc3 in the integration tests
		//               using DirichletT = typename HeatEqBundle::template DirichletBC<SourceCallableT>;
		//               auto bc = std::make_shared<fem::boundary::BoundaryCondition<DirichletT>>(
		//                   fem::boundary::BoundaryCondition<DirichletT>{
		//                       bcCfg.boundaryID,
		//                       {fem::boundary::BCCategory::Essential},
		//                       DirichletT{SourceCallableT(bcCfg.expression)}
		//                   });
		//               essentialBCs_.registerBC<DirichletT>(bc);
		//           } else if (form == config::BoundaryConditionConfig::Form::FluxBC) {
		//               // natural (flux) BC -- mirrors bc2 in tests/integration/fem/cpu/HeatEquationMinimal.cpp,
		//               // which registers with HeatEqBundle::DefaultModelBdy, not a conductivity model.
		//               using FluxFunctionT = typename HeatEqBundle::template FluxBC<FluxCallableT>;
		//               using FluxFormT = typename HeatEqBundle::template FluxForm<FluxCallableT>;
		//               fluxForms_.push_back(FluxFormsT{FluxFormT{FluxFunctionT{FluxCallableT(bcCfg.expression)}}});
		//               auto bc = std::make_shared<fem::boundary::BoundaryCondition<FluxFunctionT>>(
		//                   fem::boundary::BoundaryCondition<FluxFunctionT>{
		//                       bcCfg.boundaryID,
		//                       {fem::boundary::BCCategory::Natural},
		//                       FluxFunctionT{FluxCallableT(bcCfg.expression)}
		//                   });
		//               naturalBCs_.registerBC<FluxFunctionT>(bc, fluxForms_.back(), defaultModelBdy_);
		//           }
		//       }
		//   }

		// --- 4. constraints, sized allocation -----------------------------------------
		// MUST happen after every essential BC above is registered.
		//
		//   topoDOF_.buildConstraints<typename HeatEqBundle::Basis>(essentialBCs_);
		//
		//   if (config_.solver.linear->operatorType == config::LinearSolverConfig::OperatorType::CSR) {
		//       K_ = fem::assembly::Assembler<Backend>::template createMatrix<HeatEqBundle::NumDOFs>(mesh_, topoDOF_);
		//   }
		//   F_ = fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_);
		//   U_ = fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_);
		//   U_.zero();

		// --- 5. linear solver runner --------------------------------------------------
		// config_.solver.linear is a std::optional -- HeatDispatcher should validate
		// it's populated before ever constructing a HeatProblem; treat an empty
		// optional here as a logic error (throw), not a silently-defaulted config.
		// Matrix vs matrix-free is a runtime choice made once, right here, by picking
		// which OperatorType to hand to the factory -- see FEMOperatorT's docs in the
		// header for why the temporary Operator objects below are safe to let go out
		// of scope immediately after this call.
		//
		//   if (operatorType == CSR) {
		//       CSROperatorT op(K_);
		//       linearSolverRunner_ = solver::linear::makeLinearSolverRunner<CSROperatorT, VectorT>(
		//           op, topoDOF_.numFreeDOFs(), *config_.solver.linear, "heateq");
		//   } else {
		//       FEMOperatorT op(assembler_, mesh_, topoDOF_, /*time*/ Real(0), conductivityModel_, matrixForms_);
		//       linearSolverRunner_ = solver::linear::makeLinearSolverRunner<FEMOperatorT, VectorT>(
		//           op, topoDOF_.numFreeDOFs(), *config_.solver.linear, "heateq");
		//   }
	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::assembleSystem(Real time) {

		// Mirrors tests/integration/full/cpu/HeatEquationMinimal.cpp's
		// MatrixCGSolverBilinearSolP1 body (matrix path) / MatrixFreeCGSolver body
		// (matrix-free path skips the assembleMatrix call entirely -- K_ is never
		// touched, the operator applies DiffusionForm on the fly instead).
		//
		//   if (matrix path) {
		//       zero K_ (see Assembler.tpp/CSRMatrix for whatever the established
		//       "clear then re-fill" idiom is -- assembleMatrix accumulates with +=,
		//       it does not overwrite);
		//       fem::assembly::Assembler<Backend>::template assembleMatrix<
		//           HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol,
		//           typename HeatEqBundle::ConstantConductivityModel, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType
		//       >(mesh_, topoDOF_, time, conductivityModel_, matrixForms_, U_, K_);
		//   }
		//
		//   zero F_;
		//   fem::assembly::Assembler<Backend>::template assembleVector<
		//       HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol,
		//       typename HeatEqBundle::DefaultModel, SourceFormsT, typename HeatEqBundle::QuadratureVolumeType
		//   >(mesh_, topoDOF_, time, defaultModel_, sourceForms_, U_, F_);        // note: DefaultModel, not conductivityModel_
		//
		//   bcApplicator_.template applyNaturalBCs<
		//       HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPBdy,
		//       typename HeatEqBundle::QuadratureBoundaryType
		//   >(mesh_, topoDOF_, naturalBCs_, time, F_);
		//
		//   bcApplicator_.template applyEssentialBCs<
		//       HeatEqBundle::NumDOFs, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol,
		//       typename HeatEqBundle::ConstantConductivityModel, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType
		//   >(mesh_, topoDOF_, essentialBCs_, time, conductivityModel_, matrixForms_, F_);  // essential BC lifting needs the matrix-side model/forms, not DefaultModel
	}

	template<typename Backend, typename HeatEqBundle>
	bool HeatProblem<Backend, HeatEqBundle>::solveLinear() {

		// Mirrors the integration tests' solver.solve(report, ...) call, but through
		// linearSolverRunner_ instead of a hand-built cg::Solver -- that's the whole
		// point of the *Runner indirection (see solver::linear::LinearSolverFactory).
		//
		//   linalg::solver::SolverReport<VectorT> report;
		//   bool converged = linearSolverRunner_->solve(F_, U_, report);
		//   return converged;
		return false;
	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::writeOutput(Index step) const {

		// Mirrors tests/integration/full/cpu/HeatEquationMinimal.cpp's writeVTK call
		// at the end of both MatrixCGSolverBilinearSolP1 and MatrixFreeCGSolver.
		// config_.output gates whether/how often this actually writes anything
		// (config_.output.vtk, config_.output.writeFrequency) -- check those before
		// calling FieldIO at all.
		//
		//   if (!config_.output.vtk || (step % config_.output.writeFrequency != 0)) return;
		//
		//   const std::string filename = config_.output.directory + "/" + config_.output.prefix + "_" + std::to_string(step) + ".vtk";
		//   io::FieldIO::writeVTK<HeatEqBundle::NumDOFs>(
		//       mesh_, topoDOF_, essentialBCs_, /*time*/ Real(0), U_.data(), {"T"}, filename);
		//
		// ("T" is a placeholder field name -- the integration tests use "theta"
		// hardcoded; consider whether this should come from config instead.)
	}

	template<typename Backend, typename HeatEqBundle>
	void HeatProblem<Backend, HeatEqBundle>::writeLog() const {

		// No LoggingConfig exists anywhere in the codebase yet (grep confirms
		// HeatConfigParser never reads root["logging"], and
		// include/utils/logging/driver/ is an empty directory) -- steady.yaml's
		// logging: section is currently unconsumed. For now this can only report to
		// stdout, e.g. via the same utils::logging::ConsoleLogger pattern
		// LinearSolverFactory already uses internally. Revisit once a LoggingConfig
		// struct + parser exist.
	}

} // namespace pdesolver::application::heateq::problem
