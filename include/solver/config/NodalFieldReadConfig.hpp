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
				std::string expression;
				std::string file;
				std::string unit;

			}; // struct NodalFieldReadConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
