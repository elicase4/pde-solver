#ifndef PDESOLVER_SOLVER_CONFIG_OUTPUTCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_OUTPUTCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace pdesolver {
	namespace solver {
		namespace config {

			struct OutputConfig {

				std::string directory = "output";

				bool vtk = true;

				Index writeFrequency = 1;

				std::string prefix = "solution";

			}; // struct OutputConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
