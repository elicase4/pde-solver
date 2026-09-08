#ifndef PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP
#define PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP

#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
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
#include "fem/form/FormRegistry.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"

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

					const solver::SolverInstance& solverInstance() const { return solverInstance_; }

					const VectorT& solution() const { return *U_; }

					Index numDOFs() const { return topoDOF_.numFreeDOFs(); }
					
					Real residualNorm() const;

					void writeOutput(Index step) const;

					void writeLog() const;

					void evaluateMonitors() const;

				private:

					using SourceCallableT = utils::expression::ScalarExpression;
					using FluxCallableT = utils::expression::VectorExpression;

					using NodalScalarSourceT = typename HeatEqBundle::NodalScalarSource;
					using NodalFluxSourceT = typename HeatEqBundle::NodalFluxSource;

					using MatrixFormsT = typename HeatEqBundle::MatrixForms;
					using ExpressionSourceFormsT = typename HeatEqBundle::template ExpressionSourceForms<SourceCallableT>;
					using NodalSourceFormsT = typename HeatEqBundle::template NodalSourceForms<NodalScalarSourceT>;
					using ExpressionFluxFormsT = typename HeatEqBundle::template ExpressionFluxForms<FluxCallableT>;
					using NodalFluxFormsT = typename HeatEqBundle::template NodalFluxForms<NodalFluxSourceT>;

					using DirichletExpressionT = typename HeatEqBundle::template DirichletExpression<SourceCallableT>;
					using DirichletNodalT = typename HeatEqBundle::template DirichletNodal<NodalScalarSourceT>;
					using FluxFunctionExpressionT = typename HeatEqBundle::template FluxBCExpression<FluxCallableT>;
					using FluxFunctionNodalT = typename HeatEqBundle::template FluxBCNodal<NodalFluxSourceT>;

					template<fem::quantity::Reduction Mode>
					using MonitorQuantitiesFor = typename HeatEqBundle::template QuantityForms<typename HeatEqBundle::template ReducedQuantity<typename HeatEqBundle::HeatFluxIntegrand, Mode>>;
					using MonitorQuantitiesIntegralT = MonitorQuantitiesFor<fem::quantity::Reduction::Integral>;
					using MonitorQuantitiesAverageT = MonitorQuantitiesFor<fem::quantity::Reduction::Average>;
					using MonitorRegistryIntegralT = typename HeatEqBundle::template BoundaryQuantityRegistry<MonitorQuantitiesIntegralT>;
					using MonitorRegistryAverageT = typename HeatEqBundle::template BoundaryQuantityRegistry<MonitorQuantitiesAverageT>;
					using MonitorCombinationIntegralT = typename HeatEqBundle::template BoundaryQuantityCombination<MonitorQuantitiesIntegralT>;
					using MonitorCombinationAverageT = typename HeatEqBundle::template BoundaryQuantityCombination<MonitorQuantitiesAverageT>;

					using CSROperatorT = linalg::op::CSROperator<MatrixT>;

					template<typename ConductivityModelT>
					using FEMOperatorFor = linalg::op::FEMOperator<fem::assembly::Assembler<Backend>, topology::TopologicalDOF<HeatEqBundle::NumDOFs>, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, ConductivityModelT, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType>;

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

					typename HeatEqBundle::ConductivityModelVariant conductivityModel_;
					typename HeatEqBundle::ConductivityModelVariantBdy conductivityModelBdy_;

					typename HeatEqBundle::DefaultModel defaultModel_;
					typename HeatEqBundle::DefaultModelBdy defaultModelBdy_;

					MonitorQuantitiesIntegralT monitorQuantitiesIntegral_;
					MonitorQuantitiesAverageT monitorQuantitiesAverage_;
					mutable MonitorRegistryIntegralT monitorRegistryIntegral_;
					mutable MonitorRegistryAverageT monitorRegistryAverage_;
					std::vector<std::pair<std::string, MonitorCombinationIntegralT>> monitorCombinationsIntegral_;
					std::vector<std::pair<std::string, MonitorCombinationAverageT>> monitorCombinationsAverage_;

					MatrixFormsT matrixForms_;
					std::optional<ExpressionSourceFormsT> sourceForms_;
					std::optional<NodalSourceFormsT> nodalSourceForms_;

					std::vector<std::unique_ptr<ExpressionFluxFormsT>> expressionFluxForms_;
					std::vector<std::unique_ptr<NodalFluxFormsT>> nodalFluxForms_;

					std::unique_ptr<MatrixT> K_;
					std::unique_ptr<VectorT> F_;
					std::unique_ptr<VectorT> U_;

					std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> linearSolverRunner_;

				}; // class HeatProblem

			} // namespace problem
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#include "application/heateq/problem/HeatProblem.tpp"

#endif
