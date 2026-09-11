#ifndef PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP
#define PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP

#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "application/heateq/config/HeatConfig.hpp"

#include "core/Types.hpp"

#include "io/fieldio/FieldIO.hpp"
#include "io/fieldio/NodalFileValueSource.hpp"
#include "io/fieldio/NodalValueSourceAdapter.hpp"

#include "fem/assembly/Assembler.hpp"
#include "fem/boundary/BoundaryApplicator.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/NaturalBoundaryRegistry.hpp"
#include "fem/eval/ModelRegistry.hpp"
#include "fem/form/FormRegistry.hpp"
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

					HeatProblem(const config::HeatConfig& config, mesh::Mesh mesh, typename HeatEqBundle::Basis basis, typename HeatEqBundle::QuadratureVolumeType quadratureVolume, typename HeatEqBundle::QuadratureBoundaryType quadratureBoundary);

					HeatProblem(const HeatProblem&) = delete;
					HeatProblem& operator=(const HeatProblem&) = delete;
					HeatProblem(HeatProblem&&) = delete;
					HeatProblem& operator=(HeatProblem&&) = delete;

					void assembleSystem(Real time);

					bool solveLinear();

					bool solveLinearStep();

					void advanceTimestep();

					void setDt(Real dt);

					const solver::SolverInstance& solverInstance() const { return solverInstance_; }

					const VectorT& solution() const { return *U_; }

					Index numDOFs() const { return topoDOF_.numFreeDOFs(); }
					
					Real residualNorm() const;

					void writeOutput(Index step) const;

					void writeLog() const;

					void evaluateMonitors(Index tick, Real time) const;

				private:

					using OpType = solver::config::LinearSolverConfig::OperatorType;

					// assembly helpers
					void loadNodalField(const solver::config::NodalFieldReadConfig& cfg, VectorT& target) const;
					void assembleStiffnessMatrix(Real time);   // steady K -> K_
					void assembleTransientOperatorMatrix();     // K + M/dt -> K_
					void assembleLoad(Real time);              // volumetric source -> F_
					void addTransientMassTerm(Real time);      // F_ += (C/dt)*M*U_prev
					void applyNatural(Real time);              // F_ += flux BCs
					void applyEssential(Real time);            // constrain F_ against the active operator

					Real transientMassCoefficient() const { return Real(1) / dt_; }

					template<typename OperatorT>
					void makeLinearRunner(const OperatorT& op);

					using SourceCallableT = utils::expression::ScalarExpression;
					using FluxCallableT = utils::expression::VectorExpression;

					using NodalScalarSourceT = typename HeatEqBundle::NodalScalarSource;
					using NodalFluxSourceT = typename HeatEqBundle::NodalFluxSource;

					using StiffnessFormsT = typename HeatEqBundle::StiffnessForms;
					using ExpressionSourceFormsT = typename HeatEqBundle::template ExpressionSourceForms<SourceCallableT>;
					using NodalSourceFormsT = typename HeatEqBundle::template NodalSourceForms<NodalScalarSourceT>;
					using ExpressionFluxFormsT = typename HeatEqBundle::template ExpressionFluxForms<FluxCallableT>;
					using NodalFluxFormsT = typename HeatEqBundle::template NodalFluxForms<NodalFluxSourceT>;

					using DirichletExpressionT = typename HeatEqBundle::template DirichletExpression<SourceCallableT>;
					using DirichletNodalT = typename HeatEqBundle::template DirichletNodal<NodalScalarSourceT>;
					using FluxFunctionExpressionT = typename HeatEqBundle::template FluxBCExpression<FluxCallableT>;
					using FluxFunctionNodalT = typename HeatEqBundle::template FluxBCNodal<NodalFluxSourceT>;

					// TODO: generalize for more monitors
					template<fem::quantity::Reduction Mode>
					using MonitorQuantitiesFor = typename HeatEqBundle::template QuantityForms<typename HeatEqBundle::template ReducedQuantity<typename HeatEqBundle::HeatFluxIntegrand, Mode>>;
					using MonitorQuantitiesIntegralT = MonitorQuantitiesFor<fem::quantity::Reduction::Integral>;
					using MonitorQuantitiesAverageT = MonitorQuantitiesFor<fem::quantity::Reduction::Average>;
					using MonitorRegistryIntegralT = typename HeatEqBundle::template BoundaryQuantityRegistry<MonitorQuantitiesIntegralT>;
					using MonitorRegistryAverageT = typename HeatEqBundle::template BoundaryQuantityRegistry<MonitorQuantitiesAverageT>;
					using MonitorCombinationIntegralT = typename HeatEqBundle::template BoundaryQuantityCombination<MonitorQuantitiesIntegralT>;
					using MonitorCombinationAverageT = typename HeatEqBundle::template BoundaryQuantityCombination<MonitorQuantitiesAverageT>;

					template<typename CombinationT>
					struct MonitorOutput {
						std::string name;
						CombinationT combination;
						bool toConsole;
						std::string unit;
						utils::logging::CsvWriter csv;
					}; // struct MonitorOutput

					using CSROperatorT = linalg::op::CSROperator<MatrixT>;

					using TransientOperatorFormsT = typename HeatEqBundle::TransientOperatorForms;
					using MassFormsT = typename HeatEqBundle::MassForms;

					// matrix-free operators
					template<typename ModelT, typename FormsT>
					using MatrixFreeOperator = linalg::op::FEMOperator<fem::assembly::Assembler<Backend>, topology::TopologicalDOF<HeatEqBundle::NumDOFs>, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, ModelT, FormsT, typename HeatEqBundle::QuadratureVolumeType>;
					using SteadyMatrixFreeOperatorT = MatrixFreeOperator<typename HeatEqBundle::ConductivityModel, StiffnessFormsT>;
					using TransientMatrixFreeOperatorT = MatrixFreeOperator<typename HeatEqBundle::MaterialModel, TransientOperatorFormsT>;

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

					MonitorQuantitiesIntegralT monitorQuantitiesIntegral_;
					MonitorQuantitiesAverageT monitorQuantitiesAverage_;
					mutable MonitorRegistryIntegralT monitorRegistryIntegral_;
					mutable MonitorRegistryAverageT monitorRegistryAverage_;
					std::vector<MonitorOutput<MonitorCombinationIntegralT>> monitorOutputsIntegral_;
					std::vector<MonitorOutput<MonitorCombinationAverageT>> monitorOutputsAverage_;

					StiffnessFormsT stiffnessForms_;
					std::optional<ExpressionSourceFormsT> sourceForms_;
					std::optional<NodalSourceFormsT> nodalSourceForms_;

					std::vector<std::unique_ptr<ExpressionFluxFormsT>> expressionFluxForms_;
					std::vector<std::unique_ptr<NodalFluxFormsT>> nodalFluxForms_;

					// transient-only
					Real dt_ = Real(0);
					typename HeatEqBundle::DensityModel densityModel_;
					typename HeatEqBundle::SpecificHeatModel specificHeatModel_;
					typename HeatEqBundle::MaterialModel materialModel_;
					typename HeatEqBundle::MassMaterialModel massMaterialModel_;
					TransientOperatorFormsT transientOperatorForms_;
					MassFormsT massForms_;

					std::unique_ptr<MatrixT> K_;
					std::unique_ptr<VectorT> F_;
					std::unique_ptr<VectorT> U_;
					std::unique_ptr<VectorT> U_prev_;      // transient-only
					std::unique_ptr<VectorT> massTimesUprev_; // transient-only scratch for M*U_prev

					std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> linearSolverRunner_;

				}; // class HeatProblem

			} // namespace problem
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#include "application/heateq/problem/HeatProblem.tpp"

#endif
