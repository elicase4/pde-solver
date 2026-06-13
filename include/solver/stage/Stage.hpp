#ifndef PDESOLVER_STAGE_STAGE_HPP
#define PDESOLVER_STAGE_STAGE_HPP

#include <concepts>

namespace pdesolver {
	namespace solver {
		namespace stage {

			template<typename S>
			concept Stage = requires(S& stage) {

				{ stage.initialize() };
				{ stage.assemble() };
				{ stage.solve() } -> std::same_as<bool>;
				{ stage.finalize() };

			}; // concept Stage

		} // namespace stage
	} // namespace solver
} // namespace pdesolver

#endif
