#ifndef PDESOLVER_DRIVER_DRIVER_HPP
#define PDESOLVER_DRIVER_DRIVER_HPP

#include <concepts>

namespace pdesolver {
	namespace solver {
		namespace driver {

			template<typename D, typename StageType>
			concept Driver = requires(D& driver, StageType& stage) {
				
				{ driver.solve(stage) } -> std::same_as<bool>;

			}; // concept Driver

		} // namespace driver
	} // namespace solver
} // namespace pdesolver

#endif
