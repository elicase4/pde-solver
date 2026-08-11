#ifndef PDESOLVER_SOLVER_STAGE_NONLINEARCAPABLESTAGE_HPP
#define PDESOLVER_SOLVER_STAGE_NONLINEARCAPABLESTAGE_HPP

#include <concepts>

#include "core/Types.hpp"
#include "solver/stage/Stage.hpp"

namespace pdesolver {
	namespace solver {
		namespace stage {

			// Refines Stage with what a Newton/Picard-style nonlinear iteration needs to drive
			// its own loop without knowing anything about the underlying equation's matrix/vector
			// types. A stage that's always purely linear never needs to satisfy this -- only
			// Steady<StageT>/Transient<StageT,...> require the base Stage concept; this is
			// additionally required only where NewtonRunner/PicardRunner are instantiated.
			//
			// Contract:
			//   - assemble() (from Stage) re-evaluates the residual R(U) and tangent J(U) at the
			//     stage's current solution U.
			//   - residualNorm() reports ||R(U)|| from that last assemble() call.
			//   - solveLinearStep() solves J(U)*deltaU = -R(U) for the correction and applies
			//     U += deltaU internally (the stage owns U; nothing external is passed in or
			//     read back here). Returns whether the linear solve itself succeeded -- NOT
			//     whether the outer nonlinear iteration has converged, which the caller judges
			//     via residualNorm() against its own tolerance.
			template<typename S>
			concept NonlinearCapableStage = Stage<S> && requires(S& stage) {

				{ stage.residualNorm() } -> std::convertible_to<Real>;

				{ stage.solveLinearStep() } -> std::same_as<bool>;

			}; // concept NonlinearCapableStage

		} // namespace stage
	} // namespace solver
} // namespace pdesolver

#endif
