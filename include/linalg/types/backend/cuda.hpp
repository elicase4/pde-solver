#ifndef PDESOLVER_LINALG_TYPES_BACKEND_CUDA_HPP
#define PDESOLVER_LINALG_TYPES_BACKEND_CUDA_HPP

#include <cuda_runtime.h>
#include <memory>
#include <stdexcept>

#include "core/Types.hpp"

namespace pdesolver {
	namespace linalg {
		namespace types {
			namespace backend {
				
				template<typename T>
				struct CudaDeleter {
					void operator()(T* ptr) const noexcept {
						if (ptr) cudaFree(ptr);
					}
				}; // struct CudaDeleter

				class CUDA {
				public:
					
					enum class Direction {HostToDevice, DeviceToHost, DeviceToDevice};
					
					template<typename T>
					using Ptr = std::unique_ptr<T, CudaDeleter<T>>;

					template<typename T>
					static Ptr<T> alloc(Index n){
						T* ptr = nullptr;
						cudaError_t err = cudaMalloc(&ptr, n*sizeof(T));
						if (err != cudaCuccess){
							throw std::runtime_error("cudaMalloc failure");
						}
						cudaMemset(ptr, 0, n*sizeof(T));
						return Ptr<T>(ptr);
					}

					template<typename T>
					static void copy(T* dst, const T* src, Index n, Direction dir){
						cudaMemcpyKind kind = cudaMemcpyDefault;
						switch(dir) {
							case (Direction::HostToDevice){
								kind = cudaMemcpyHostToDevice;
								break;
							}
							case (Direction::DeviceToHost){
								kind = cudaMemcpyDeviceToHost;
								break;
							}
							case (Direction::DeviceToDevice){
								kind = cudaMemcpyDeviceToDevice;
								break;
							}
						}
						cudaMemcpy(dst, src, n*sizeof(T), kind);
					}
				
				}; // class CUDA

			} // namespace backend
		} // namespace types
	} // namespace linalg
} // namespace pdesolver

#endif
