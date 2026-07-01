#ifndef PDESOLVER_APPLICATION_HEATEQ_BOUNDARYCONDITIONCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_BOUNDARYCONDITIONCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			struct BoundaryConditionConfig {

				enum class Type {
					Value,
					Flux
				}; // enum class Type

				Int boundaryID;

				Type type;

				std::string expression;

			}; // struct BoundaryConditionConfig

		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
