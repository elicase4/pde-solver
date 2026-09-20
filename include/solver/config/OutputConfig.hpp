#ifndef PDESOLVER_SOLVER_CONFIG_OUTPUTCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_OUTPUTCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace pdesolver {
	namespace solver {
		namespace config {

			struct OutputConfig {

				enum class Format { VTK, VTU };

				std::string directory = "output";

				Format format = Format::VTK;

				Index writeFrequency = 1;

				std::string prefix = "solution";

			}; // struct OutputConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
