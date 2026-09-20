#ifndef PDESOLVER_SOLVER_PROBLEM_PROBLEM_HPP
#define PDESOLVER_SOLVER_PROBLEM_PROBLEM_HPP

#include <concepts>
#include <memory>

#include "core/Types.hpp"

#include "fem/assembly/ElementMap.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"

namespace pdesolver {
	namespace solver {
		namespace problem {

			template<typename P>
			concept Problem = requires(P& p, Real time, typename P::MatrixT& K, typename P::VectorT& V) {

				typename P::VectorT;
				typename P::MatrixT;

				{ P::NumDOFs } -> std::convertible_to<Index>;
				{ P::NumAuxStates } -> std::convertible_to<Index>;

				p.stiffnessForms();
				p.stiffnessModel();
				p.tangentStiffnessForms();

				{ p.createMatrix() } -> std::same_as<typename P::MatrixT>;
				{ p.createVector() } -> std::same_as<typename P::VectorT>;

				{ p.assembleLoad(time) };
				{ p.applyNatural(time) };

				{ p.template assembleMatrix<fem::assembly::GatherMode::Free>(time, p.stiffnessForms(), p.stiffnessModel(), {}, K) };
				{ p.template assembleVector<fem::assembly::GatherMode::Free>(time, p.stiffnessForms(), p.stiffnessModel(), V, {}, V) };
				{ p.template assembleResidual<fem::assembly::GatherMode::Free>(time, p.stiffnessForms(), p.stiffnessModel(), {}, V) };
				{ p.applyEssential(time, p.stiffnessForms(), p.stiffnessModel(), V) };

				{ p.template makeLinearRunner<fem::assembly::GatherMode::Free>(&time, p.stiffnessForms(), p.stiffnessModel(), nullptr, {}, K) } -> std::same_as<std::unique_ptr<linalg::solver::LinearSolverRunner<typename P::VectorT>>>;

				{ p.numFreeDOFs() } -> std::convertible_to<Index>;
				{ p.solverInstance() };
				{ p.equationLabel() };
				{ p.loggingConfig() };

				{ p.U() } -> std::convertible_to<typename P::VectorT&>;
				{ p.F() } -> std::convertible_to<typename P::VectorT&>;

			}; // concept Problem

		} // namespace problem
	} // namespace solver
} // namespace pdesolver

#endif
