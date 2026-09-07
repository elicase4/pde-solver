#ifndef PDESOLVER_APPLICATION_HEATEQ_PARSER_SOURCECONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_PARSER_SOURCECONFIG_HPP

#include "solver/config/NodalFieldReadConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace config {

				struct SourceConfig {

					enum class Type {
						VolumetricHeatSource
					}; // enum class Type

					Type type;

					// the heat equation has one DOF, so NodalFieldReadConfig's scalar
					// expression shape fits directly, same as InitialConditionConfig
					solver::config::NodalFieldReadConfig read;

				}; // struct SourceConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
