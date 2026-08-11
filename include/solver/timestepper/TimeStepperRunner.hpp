#ifndef PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPERRUNNER_HPP
#define PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPERRUNNER_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace solver {
		namespace timestepper {

			// Type-erased time integration scheme, mirroring LinearSolverRunner/
			// NonlinearSolverRunner's role: TimeStepperConfig::Type (ForwardEuler/BackwardEuler/
			// GeneralizedAlpha/RK4) is a runtime YAML value, not a hot-loop-frequency choice (it's
			// selected once per timestep, nowhere near quadrature-point frequency), so it should
			// be resolved via virtual dispatch the same way LinearSolverConfig::Type already is --
			// not folded into fem::dispatch's compile-time combinatorics.
			//
			// A concrete implementation (e.g. a future BackwardEulerRunner) composes -- does not
			// parallel -- LinearSolverRunner/NonlinearSolverRunner: it holds a reference to
			// whichever one applies for the wrapped stage's resolved SolverMode, and calls into it
			// once per step to solve that step's system.
			//
			// TODO: advance()'s (U, U_prev) parameters predate the NonlinearCapableStage design
			// (where the stage owns U internally); revisit whether these should drop in favor of
			// the stage-owns-everything pattern once a concrete stepper is actually implemented.
			template<typename VectorType>
			class TimeStepperRunner {
			public:

				virtual ~TimeStepperRunner() = default;

				virtual bool advance(VectorType& U, VectorType& U_prev) = 0;

				virtual bool finished() const = 0;

				virtual Real currentTime() const = 0;

			}; // class TimeStepperRunner

		} // namespace timestepper
	} // namespace solver
} // namespace pdesolver

#endif
