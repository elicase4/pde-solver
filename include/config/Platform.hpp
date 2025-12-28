#ifndef PDESOLVER_PLATFORM_HPP
#define PDESOLVER_PLATFORM_HPP

#ifdef __CUDACC__
	#define PDE_HOST __host__
	#define PDE_DEVICE __device__
	#define PDE_INLINE __forceinline__
#else
	#define PDE_HOST
	#define PDE_DEVICE
	#define PDE_INLINE inline
#endif

#endif
