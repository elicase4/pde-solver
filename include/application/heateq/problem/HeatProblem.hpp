#ifndef PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP
#define PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "application/heateq/config/HeatConfig.hpp"

#include "core/Types.hpp"

#include "io/visualization/VisualizationWriter.hpp"
#include "io/fieldio/NodalFileValueSource.hpp"
#include "io/fieldio/NodalValueSourceAdapter.hpp"

#include "fem/assembly/Assembler.hpp"
#include "fem/boundary/BoundaryApplicator.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/NaturalBoundaryRegistry.hpp"
#include "fem/form/FormRegistry.hpp"
#include "fem/quantity/MonitorGroup.hpp"
#include "fem/quantity/MonitorGroupImpl.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"
#include "fem/quantity/QuantityUnits.hpp"

#include "utils/logging/core/CsvWriter.hpp"

#include "linalg/operations/VectorOps.hpp"
#include "linalg/operator/CSROperator.hpp"
#include "linalg/operator/FEMOperator.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"

#include "mesh/Mesh.hpp"

#include "solver/SolverInstance.hpp"
#include "solver/config/NodalFieldReadConfig.hpp"
#include "solver/linear/LinearSolverFactory.hpp"
#include "solver/logging/LoggerFactory.hpp"

#include "topology/TopologicalDOF.hpp"

#include "utils/expression/ScalarExpression.hpp"
#include "utils/expression/VectorExpression.hpp"
#include "utils/logging/driver/Logger.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace problem {

				template<typename Backend, typename HeatEqBundle>
				class HeatProblem {
				public:

					using VectorT = linalg::types::Vector<Real, Backend>;
					using MatrixT = linalg::types::CSRMatrix<Real, Backend>;

					static constexpr Index NumDOFs = HeatEqBundle::NumDOFs;
					static constexpr Index NumAuxStates = HeatEqBundle::EvalQPVol::NumAuxStates;

					HeatProblem(const config::HeatConfig& config, mesh::Mesh mesh, typename HeatEqBundle::Basis basis, typename HeatEqBundle::QuadratureVolumeType quadratureVolume, typename HeatEqBundle::QuadratureBoundaryType quadratureBoundary);

					HeatProblem(const HeatProblem&) = delete;
					HeatProblem& operator=(const HeatProblem&) = delete;
					HeatProblem(HeatProblem&&) = delete;
					HeatProblem& operator=(HeatProblem&&) = delete;

					// stateless
					typename HeatEqBundle::StiffnessForms stiffnessForms() const { return {}; }
					typename HeatEqBundle::MassForms massForms() const { return {}; }
					typename HeatEqBundle::TangentDiffusionForms tangentStiffnessForms() const { return {}; }
					typename HeatEqBundle::TangentMassForms tangentMassForms() const { return {}; }

					const typename HeatEqBundle::ConductivityModel& stiffnessModel() const { return conductivityModel_; }
					const typename HeatEqBundle::MassModel& massModel() const { return massModel_; }

					MatrixT createMatrix() const { return fem::assembly::Assembler<Backend>::template createMatrix<HeatEqBundle::NumDOFs>(mesh_, topoDOF_); }
					VectorT createVector() const { return fem::assembly::Assembler<Backend>::template createVector<HeatEqBundle::NumDOFs>(mesh_, topoDOF_); }

					void assembleLoad(Real time);
					void applyNatural(Real time);

					template<fem::assembly::GatherMode Mode, typename FormsT, typename ModelT>
					void assembleMatrix(Real time, const FormsT& forms, const ModelT& model, const std::array<const VectorT*, NumAuxStates>& auxStates, MatrixT& K);

					template<fem::assembly::GatherMode Mode, typename FormsT, typename ModelT>
					void assembleVector(Real time, const FormsT& forms, const ModelT& model, const VectorT& gatherSource, const std::array<const VectorT*, NumAuxStates>& auxStates, VectorT& V);

					template<fem::assembly::GatherMode Mode, typename FormsT, typename ModelT>
					void assembleResidual(Real time, const FormsT& forms, const ModelT& model, const std::array<const VectorT*, NumAuxStates>& auxStates, VectorT& R);

					template<typename FormsT, typename ModelT>
					void applyEssential(Real time, const FormsT& forms, const ModelT& model, VectorT& F);

					template<fem::assembly::GatherMode Mode, typename FormsT, typename ModelT>
					std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> makeLinearRunner(const Real* time, const FormsT& forms, const ModelT& model, const VectorT* fieldSource, const std::array<const VectorT*, NumAuxStates>& auxStates, MatrixT& K);

					Index numFreeDOFs() const { return topoDOF_.numFreeDOFs(); }
					const solver::SolverInstance& solverInstance() const { return solverInstance_; }
					std::string equationLabel() const { return "Heat Equation"; }
					const solver::config::LoggingConfig& loggingConfig() const { return config_.logging; }

					VectorT& U() { return *U_; }
					VectorT& F() { return *F_; }
					VectorT& U_prev() { return *U_prev_; }

					const VectorT& solution() const { return *U_; }

					void writeOutput(Index step, Real time) const;

					void writeLog() const;

					void evaluateMonitors(Index tick, Real time) const;

				private:

					using OpType = solver::config::LinearSolverConfig::OperatorType;

					void loadNodalField(const solver::config::NodalFieldReadConfig& cfg, VectorT& target) const;

					using ExpressionSourceFormsT = typename HeatEqBundle::template ExpressionSourceForms<utils::expression::ScalarExpression>;
					using NodalSourceFormsT = typename HeatEqBundle::template NodalSourceForms<typename HeatEqBundle::NodalScalarSource>;
					using ExpressionFluxFormsT = typename HeatEqBundle::template ExpressionFluxForms<utils::expression::VectorExpression>;
					using NodalFluxFormsT = typename HeatEqBundle::template NodalFluxForms<typename HeatEqBundle::NodalFluxSource>;

					using DirichletExpressionT = typename HeatEqBundle::template DirichletExpression<utils::expression::ScalarExpression>;
					using DirichletNodalT = typename HeatEqBundle::template DirichletNodal<typename HeatEqBundle::NodalScalarSource>;
					using FluxFunctionExpressionT = typename HeatEqBundle::template FluxBCExpression<utils::expression::VectorExpression>;
					using FluxFunctionNodalT = typename HeatEqBundle::template FluxBCNodal<typename HeatEqBundle::NodalFluxSource>;

					struct MonitorOutput {
						std::string name;
						bool toConsole;
						std::string unit;
						utils::logging::CsvWriter csv;
					}; // struct MonitorOutput

					std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> makeMonitorGroup(config::MonitorConfig::Quantity quantity, fem::quantity::Reduction mode) const;

					template<typename Form>
					std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> makeMonitorGroupForForm(fem::quantity::Reduction mode) const;

					template<typename Form, fem::quantity::Reduction Mode>
					std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> makeMonitorGroupFor() const;

					config::HeatConfig config_;

					solver::SolverInstance solverInstance_;

					mesh::Mesh mesh_;
					topology::TopologicalDOF<HeatEqBundle::NumDOFs> topoDOF_;

					utils::logging::driver::Logger driverLogger_;

					typename HeatEqBundle::EvalEle evalEleTemplate_;
					typename HeatEqBundle::QuadratureVolumeType quadratureVolume_;
					typename HeatEqBundle::QuadratureBoundaryType quadratureBoundary_;

					fem::assembly::Assembler<Backend> assembler_;

					fem::boundary::BoundaryApplicator<Backend> bcApplicator_;

					fem::boundary::EssentialBoundaryRegistry essentialBCs_;
					fem::boundary::NaturalBoundaryRegistry<typename HeatEqBundle::EvalQPBdy> naturalBCs_;

					typename HeatEqBundle::ConductivityModel conductivityModel_;
					typename HeatEqBundle::ConductivityModelBdy conductivityModelBdy_;

					typename HeatEqBundle::DefaultModel defaultModel_;
					typename HeatEqBundle::DefaultModelBdy defaultModelBdy_;

					mutable std::vector<std::unique_ptr<fem::quantity::MonitorGroup>> monitorGroups_;
					std::vector<std::vector<MonitorOutput>> monitorOutputsByGroup_;
					std::vector<std::string> groupUnits_;

					std::optional<ExpressionSourceFormsT> sourceForms_;
					std::optional<NodalSourceFormsT> nodalSourceForms_;

					std::vector<std::unique_ptr<ExpressionFluxFormsT>> expressionFluxForms_;
					std::vector<std::unique_ptr<NodalFluxFormsT>> nodalFluxForms_;

					typename HeatEqBundle::DensityModel densityModel_;
					typename HeatEqBundle::SpecificHeatModel specificHeatModel_;
					typename HeatEqBundle::MassModel massModel_;

					std::unique_ptr<VectorT> F_;
					std::unique_ptr<VectorT> U_;
					std::unique_ptr<VectorT> U_prev_;

					mutable std::optional<io::visualization::VisualizationWriter> vizWriter_;

				}; // class HeatProblem

			} // namespace problem
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#include "application/heateq/problem/HeatProblem.tpp"

#endif
