#ifndef PDESOLVER_SOLVER_STAGE_TRANSIENTCAPABLESTAGE_HPP
#define PDESOLVER_SOLVER_STAGE_TRANSIENTCAPABLESTAGE_HPP

#include <concepts>

#include "core/Types.hpp"
#include "solver/stage/Stage.hpp"

namespace pdesolver {
	namespace solver {
		namespace stage {

			// Refines Stage with what a timestepper needs to drive a transient solve. Only
			// Transient<StageT>/the concrete TimeStepperRunners require this; a steady-only stage
			// never needs to satisfy it.
			//
			// Contract:
			//   - setDt(dt): the step size about to be taken. Rebuilds the effective operator
			//     (K + M/dt) if dt actually changed since the last call -- a no-op otherwise, so
			//     a Constant-step-size policy costs nothing beyond the first call.
			//   - setTime(t): the instant assemble() (from Stage) should evaluate at -- for
			//     Backward Euler this is t^{n+1}.
			//   - assemble() then builds that step's system (e.g. F^{n+1} + (C/dt) M U^n, with the
			//     effective operator K + M/dt already in place).
			//   - solve() (from Stage) solves it, leaving the result in the stage's U.
			//   - advance(): roll the step-to-step state forward (U_prev <- U).
			//   - onStepComplete(step, time): per-step output side effects (VTK, monitors).
			template<typename S>
			concept TransientCapableStage = Stage<S> && requires(S& stage, Real t, Index step) {

				{ stage.setDt(t) };
				{ stage.setTime(t) };
				{ stage.advance() };
				{ stage.onStepComplete(step, t) };

			}; // concept TransientCapableStage

		} // namespace stage
	} // namespace solver
} // namespace pdesolver

#endif
