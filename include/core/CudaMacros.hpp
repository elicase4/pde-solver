#ifndef PDESOLVER_CUDAMACROS_HPP
#define PDESOLVER_CUDAMACROS_HPP

#ifdef __CUDACC__
	#define HOST_DEVICE __host__ __device__
#else
	#define HOST_DEVICE
#endif

#ifdef __CUDACC__
	#define HOST_DEVICE_INLINE __host__ __device__ __forceinline__
#else
	#define HOST_DEVICE_INLINE inline
#endif

#endif
