#ifndef PDESOLVER_ASSEMBLER_BACKEND_TYPES
#define PDESOLVER_ASSEMBLER_BACKEND_TYPES

namespace pdesolver {
	namespace fem {
		namespace assembly {
			namespace backend {

				struct CPU {};
				struct CUDA {};
				struct MPI {};

			} // namespace backend
		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#endif
