#ifndef PDESOLVER_SOLVER_CONFIG_DISCRETIZATIONCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_DISCRETIZATIONCONFIG_HPP

#include "core/Types.hpp"
#include "fem/dof/DOFOrdering.hpp"

namespace pdesolver {
	namespace solver {
		namespace config {

			struct DiscretizationConfig {

				struct Quadrature {

					Index xi = 2;

					Index eta = 2;

					Index zeta = 2;

				}; // struct Quadrature

				Quadrature quadrature;

				fem::dof::DOFOrdering dofOrdering = fem::dof::DOFOrdering::Interleaved;

			}; // struct DiscretizationConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif

