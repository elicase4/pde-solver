#ifndef PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPER_HPP
#define PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPER_HPP

#include <concepts>

namespace pdesolver {
	namespace solver {

		template<typename T, typename VectorType>
		concept TimeStepper = requires(T& ts, Real t, Real dt, VectorType& U, const VectorType& U_prev) {
			{ ts.step(t, dt, U, U_prev) } -> std::same_as<void>;
		}; // concept TimeStepper

	} // namespace solver
} // namespace pdesolver

#endif
