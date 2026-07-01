#ifndef PDESOLVER_APPLICATION_HEATEQ_CONDUCTIVITYCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_CONDUCTIVITYCONFIG_HPP

#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			struct ConductivityConfig {

				enum class Type {
					Constant,
					Anisotropic
				}; // enum class Type

				Type type;

				Real value;

				std::vector<Real> tensor;

			}; // struct ConductivityConfig

		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
