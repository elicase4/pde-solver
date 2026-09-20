#ifndef PDESOLVER_SOLVER_STAGE_TRANSIENTCAPABLESTAGE_HPP
#define PDESOLVER_SOLVER_STAGE_TRANSIENTCAPABLESTAGE_HPP

#include <concepts>

#include "core/Types.hpp"
#include "solver/stage/Stage.hpp"

namespace pdesolver {
	namespace solver {
		namespace stage {

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
