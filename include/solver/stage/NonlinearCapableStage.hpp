#ifndef PDESOLVER_SOLVER_STAGE_NONLINEARCAPABLESTAGE_HPP
#define PDESOLVER_SOLVER_STAGE_NONLINEARCAPABLESTAGE_HPP

#include <concepts>

#include "core/Types.hpp"
#include "solver/stage/Stage.hpp"

namespace pdesolver {
	namespace solver {
		namespace stage {

			template<typename S>
			concept NonlinearCapableStage = Stage<S> && requires(S& stage) {

				{ stage.residualNorm() } -> std::convertible_to<Real>;

				{ stage.solveLinearStep() } -> std::same_as<bool>;

			}; // concept NonlinearCapableStage

		} // namespace stage
	} // namespace solver
} // namespace pdesolver

#endif
