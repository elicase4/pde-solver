#ifndef PDESOLVER_SOLVER_CONFIG_DISCRETIZATIONCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_DISCRETIZATIONCONFIG_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace solver {
		namespace config {

			struct DiscretizationConfig {

				Index quadraturePointXi = 2;

				Index quadraturePointEta = 2;

				Index quadraturePointZeta = 2;

				bool blockDOFOrdering = true;

			}; // struct DiscretizationConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif

