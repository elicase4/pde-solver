#ifndef PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPER_HPP
#define PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPER_HPP

#include <concepts>

#include "core/Types.hpp"

namespace pdesolver {
	namespace solver {

		template<typename T>
		concept TimeStepper = requires(T& ts) {

			{ ts.step() } -> std::same_as<bool>;

			{ ts.finished() } -> std::same_as<bool>;

			{ ts.currentStep() } -> std::same_as<Index>;

			{ ts.currentTime() } -> std::same_as<Real>;

		}; // concept TimeStepper

	} // namespace solver
} // namespace pdesolver

#endif
