#ifndef PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPER_HPP
#define PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPER_HPP

#include <concepts>

#include "core/Types.hpp"

namespace pdesolver {
	namespace solver {

		template<typename T, typename VectorType>
		concept TimeStepper = requires(T& ts, VectorType& U, const VectorType& U_prev) {
			
			{ ts.finished() } -> std::same_as<bool>;

			{ ts.advance(U, U_prev) } -> std::same_as<void>;

			{ ts.time() } -> std::same_as<Real>;
		
			{ ts.dt() } -> std::same_as<Real>;
		
		}; // concept TimeStepper

	} // namespace solver
} // namespace pdesolver

#endif
