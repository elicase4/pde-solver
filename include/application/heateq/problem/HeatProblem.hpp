#ifndef PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP
#define PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP

#include <memory>
#include <type_traits>
#include <variant>
#include <vector>

#include "application/heateq/config/HeatConfig.hpp"

#include "core/Types.hpp"

#include "io/FieldIO.hpp"

#include "fem/assembly/Assembler.hpp"
#include "fem/boundary/BoundaryApplicator.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/NaturalBoundaryRegistry.hpp"
#include "fem/form/FormRegistry.hpp"

#include "linalg/operator/CSROperator.hpp"
#include "linalg/operator/FEMOperator.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"

#include "mesh/Mesh.hpp"

#include "solver/SolverInstance.hpp"
#include "solver/linear/LinearSolverFactory.hpp"

#include "topology/TopologicalDOF.hpp"

#include "utils/expression/ScalarExpression.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace problem {

				template<typename Backend, typename HeatEqBundle>
				class HeatProblem {
				public:

					using VectorT = linalg::types::Vector<Real, Backend>;
					using MatrixT = linalg::types::CSRMatrix<Real, Backend>;

					// basis/quadratureVolume/quadratureBoundary are the runtime instances
					// resolved once by fem::dispatch::dispatch() from the mesh's basis
					// order and the config's quadrature-point counts (see the
					// runtime-dispatch refactor) -- HeatProblem holds them for the
					// lifetime of the problem and forwards them, unchanged, into every
					// Assembler/BoundaryApplicator/FEMOperator call.
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

				private:

					using SourceCallableT = utils::expression::ScalarExpression;
					using FluxCallableT = utils::expression::ScalarExpression;

					using MatrixFormsT = fem::form::FormRegistry<typename HeatEqBundle::DiffusionForm>;
					using SourceFormsT = fem::form::FormRegistry<typename HeatEqBundle::template SourceForm<SourceCallableT>>;
					using FluxFormsT = fem::form::FormRegistry<typename HeatEqBundle::template FluxForm<FluxCallableT>>;

					using CSROperatorT = linalg::op::CSROperator<MatrixT>;

					// Templated on the conductivity model, not fixed to ConstantConductivityModel:
					// the concrete FEMOperator type depends on which alternative of
					// HeatEqBundle::ConductivityModelVariant is active, resolved once per
					// construction/assembly via std::visit in the .tpp -- not a fem::dispatch axis,
					// since this stays entirely inside one HeatEqBundle instantiation (see the
					// conductivity-model extensibility discussion: this is what keeps the ~2000
					// compile-time dispatch instantiations from multiplying per model).
					template<typename ConductivityModelT>
					using FEMOperatorFor = linalg::op::FEMOperator<fem::assembly::Assembler<Backend>, topology::TopologicalDOF<HeatEqBundle::NumDOFs>, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, ConductivityModelT, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType>;

					config::HeatConfig config_;

					solver::SolverInstance solverInstance_;

					mesh::Mesh mesh_;
					topology::TopologicalDOF<HeatEqBundle::NumDOFs> topoDOF_;

					// evalEleTemplate_ holds the resolved Basis (order fixed for this
					// problem's lifetime); copied and rebound per-element inside
					// Assembler/BoundaryApplicator, never mutated here.
					typename HeatEqBundle::EvalEle evalEleTemplate_;
					typename HeatEqBundle::QuadratureVolumeType quadratureVolume_;
					typename HeatEqBundle::QuadratureBoundaryType quadratureBoundary_;

					fem::assembly::Assembler<Backend> assembler_;

					fem::boundary::BoundaryApplicator<Backend> bcApplicator_;
					
					fem::boundary::EssentialBoundaryRegistry essentialBCs_;
					fem::boundary::NaturalBoundaryRegistry<typename HeatEqBundle::EvalQPBdy> naturalBCs_;

					typename HeatEqBundle::ConductivityModelVariant conductivityModel_;

					typename HeatEqBundle::DefaultModel defaultModel_;
					typename HeatEqBundle::DefaultModelBdy defaultModelBdy_;

					MatrixFormsT matrixForms_;
					SourceFormsT sourceForms_;

					std::vector<std::unique_ptr<FluxFormsT>> fluxForms_;

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
