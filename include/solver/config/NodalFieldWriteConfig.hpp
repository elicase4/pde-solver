#ifndef PDESOLVER_SOLVER_CONFIG_NODALFIELDWRITECONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_NODALFIELDWRITECONFIG_HPP

#include <string>

namespace pdesolver {
	namespace solver {
		namespace config {

			// optional one-shot export of a resolved value -- unlike NodalFieldReadConfig, a
			// whole-mesh-shaped export doesn't depend on the source's expression shape, so this
			// is genuinely shared across IC/BC/source without a scalar-vs-vector split.
			struct NodalFieldWriteConfig {

				std::string file; // empty = write disabled

			}; // struct NodalFieldWriteConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
