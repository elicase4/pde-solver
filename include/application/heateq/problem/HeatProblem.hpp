#ifndef PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP
#define PDESOLVER_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP

#include <memory>
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

					HeatProblem(const config::HeatConfig& config, mesh::Mesh mesh);

					HeatProblem(const HeatProblem&) = delete;
					HeatProblem& operator=(const HeatProblem&) = delete;
					HeatProblem(HeatProblem&&) = delete;
					HeatProblem& operator=(HeatProblem&&) = delete;

					void assembleSystem(Real time);

					bool solveLinear();

					void writeOutput(Index step) const;

					void writeLog() const;

					const VectorT& solution() const { return *U_; }

					Index numDOFs() const { return topoDOF_.numFreeDOFs(); }

				private:

					using SourceCallableT = utils::expression::ScalarExpression;
					using FluxCallableT = utils::expression::ScalarExpression;

					using MatrixFormsT = fem::form::FormRegistry<typename HeatEqBundle::DiffusionForm>;
					using SourceFormsT = fem::form::FormRegistry<typename HeatEqBundle::template SourceForm<SourceCallableT>>;
					using FluxFormsT = fem::form::FormRegistry<typename HeatEqBundle::template FluxForm<FluxCallableT>>;

					using CSROperatorT = linalg::op::CSROperator<MatrixT>;
					using FEMOperatorT = linalg::op::FEMOperator<fem::assembly::Assembler<Backend>, topology::TopologicalDOF<HeatEqBundle::NumDOFs>, typename HeatEqBundle::EvalEle, typename HeatEqBundle::EvalQPVol, typename HeatEqBundle::ConstantConductivityModel, MatrixFormsT, typename HeatEqBundle::QuadratureVolumeType>;

					config::HeatConfig config_;

					mesh::Mesh mesh_;
					topology::TopologicalDOF<HeatEqBundle::NumDOFs> topoDOF_;

					fem::assembly::Assembler<Backend> assembler_;

					fem::boundary::BoundaryApplicator<Backend> bcApplicator_;
					fem::boundary::EssentialBoundaryRegistry essentialBCs_;
					fem::boundary::NaturalBoundaryRegistry<typename HeatEqBundle::EvalQPBdy> naturalBCs_;

					typename HeatEqBundle::ConstantConductivityModel conductivityModel_;

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
