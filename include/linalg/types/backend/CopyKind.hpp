#ifndef PDESOLVER_LINALG_TYPES_BACKEND_COPYKIND_HPP
#define PDESOLVER_LINALG_TYPES_BACKEND_COPYKIND_HPP

namespace pdesolver {
	namespace linalg {
		namespace types {
			namespace backend {
				
				enum class CopyKind {
					HostToHost,
					HostToDevice,
					DeviceToHost,
					DeviceToDevice,
				}; // class CopyKind
				
			} // namespace backend
		} // namespace types
	} // namespace linalg
} // namespace pdesolver

#endif
