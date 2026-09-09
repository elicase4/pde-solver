#ifndef PDESOLVER_SOLVER_CONFIG_NODALFIELDREADCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_NODALFIELDREADCONFIG_HPP

#include <string>

namespace pdesolver {
	namespace solver {
		namespace config {

			struct NodalFieldReadConfig {

				enum class Mode { 
					Expression, 
					File 
				}; // enum class Mode

				Mode mode = Mode::Expression;
				std::string expression; // only used when mode == Expression
				std::string file;       // only used when mode == File
				std::string unit;       // required, SI unit label -- e.g. "K", "W/m^3"

			}; // struct NodalFieldReadConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
