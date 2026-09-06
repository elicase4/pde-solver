#ifndef PDESOLVER_APPLICATION_HEATEQ_CONFIG_INITIALCONDITIONCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_CONFIG_INITIALCONDITIONCONFIG_HPP

#include "solver/config/NodalFieldReadConfig.hpp"
#include "solver/config/NodalFieldWriteConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace config {

				struct InitialConditionConfig {

					solver::config::NodalFieldReadConfig read;
					solver::config::NodalFieldWriteConfig write;

				}; // struct InitialConditionConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
