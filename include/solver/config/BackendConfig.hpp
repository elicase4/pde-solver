#ifndef PDESOLVER_SOLVER_CONFIG_BACKENDCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_BACKENDCONFIG_HPP

namespace pdesolver {
	namespace solver {
		namespace config {

			struct BackendConfig {

				enum class Type {
					CPU,
					CUDA
				}; // enum class Type

				Type type = Type::CPU;

			}; // struct BackendConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
