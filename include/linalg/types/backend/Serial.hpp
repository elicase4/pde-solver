#ifndef PDESOLVER_LINALG_TYPES_BACKEND_SERIAL_HPP
#define PDESOLVER_LINALG_TYPES_BACKEND_SERIAL_HPP

#include <memory>
#include <algorithm>

#include "core/Types.hpp"
#include "CopyKind.hpp"

namespace pdesolver {
	namespace linalg {
		namespace types {
			
			namespace backend {
				
				class Serial {
				public:
					template<typename T>
					using Ptr = std::unique_ptr<T[]>;

					template<typename T>
					static Ptr<T> alloc(Index n){
						return std::make_unique<T[]>(n);
					}

					template<typename T>
					static void copy(T* dst, const T* src, Index n, CopyKind){
						std::copy(src, src + n, dst);
					}

				}; // class Serial

			} // namespace backend
		} // namespace types
	} // namespace linalg
} // namespace pdesolver

#endif
