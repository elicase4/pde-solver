#ifndef PDESOLVER_APPLICATION_HEATEQ_CONFIG_MONITORCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_CONFIG_MONITORCONFIG_HPP

#include <string>
#include <vector>

#include "core/Types.hpp"

#include "fem/quantity/Reduction.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace config {

				struct MonitorTermConfig {
					Int boundary;
					Real coefficient = 1.0;
				}; // struct MonitorTermConfig

				struct MonitorOutputConfig {
					bool console = true;
					std::string file; // optional CSV path
				}; // struct MonitorOutputConfig

				struct MonitorConfig {

					enum class Quantity {
						HeatFlux
					}; // enum class Quantity

					std::string name;
					Quantity quantity;

					fem::quantity::Reduction reduction;

					std::vector<MonitorTermConfig> terms;

					MonitorOutputConfig output;

				}; // struct MonitorConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
