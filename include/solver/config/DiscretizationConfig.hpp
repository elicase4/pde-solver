#ifndef PDESOLVER_SOLVER_CONFIG_DISCRETIZATIONCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_DISCRETIZATIONCONFIG_HPP

#include "core/Types.hpp"
#include "fem/dof/DOFOrdering.hpp"

namespace pdesolver {
	namespace solver {
		namespace config {

			struct DiscretizationConfig {

				enum class BasisType {
					Lagrange
				}; // enum class BasisType

				struct Basis {

					BasisType type = BasisType::Lagrange;

					// NOTE: px/py/pz assume a single scalar field (matches the current heateq
					// use case). pz is unused on a 2D quad mesh and only takes effect for a hex
					// mesh. Revisit when a multi-field basis (e.g. Taylor-Hood for RANS) needs
					// px/py/pz to vary per field.
					Index px = 1;

					Index py = 1;

					Index pz = 1;

				}; // struct Basis

				struct Quadrature {

					Index xi = 2;

					Index eta = 2;

					Index zeta = 2;

				}; // struct Quadrature

				Basis basis;

				Quadrature quadrature;

				fem::dof::DOFOrdering dofOrdering = fem::dof::DOFOrdering::Interleaved;

			}; // struct DiscretizationConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif

