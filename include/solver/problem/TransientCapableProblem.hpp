#ifndef PDESOLVER_SOLVER_PROBLEM_TRANSIENTCAPABLEPROBLEM_HPP
#define PDESOLVER_SOLVER_PROBLEM_TRANSIENTCAPABLEPROBLEM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "solver/problem/Problem.hpp"

namespace pdesolver {
	namespace solver {
		namespace problem {

			template<typename P>
			concept TransientCapableProblem = Problem<P> && requires(P& p) {

				p.massForms();
				p.massModel();
				p.tangentMassForms();

				{ p.U_prev() } -> std::convertible_to<typename P::VectorT&>;

			}; // concept TransientCapableProblem

		} // namespace problem
	} // namespace solver
} // namespace pdesolver

#endif
