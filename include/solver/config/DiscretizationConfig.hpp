#ifndef PDESOLVER_SOLVER_CONFIG_DISCRETIZATIONCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_DISCRETIZATIONCONFIG_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace solver {
		namespace config {

			struct DiscretizationConfig {

				Index quadraturePointXi;

				Index quadraturePointEta;

				Index quadraturePointZeta;

				bool blockDOFOrdering;

			}; // struct DiscretizationConfig
		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif

