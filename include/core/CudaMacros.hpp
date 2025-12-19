#ifndef PDESOLVER_CUDAMACROS_HPP
#define PDESOLVER_CUDAMACROS_HPP

#ifdef __CUDACC__
	#define HOST_DEVICE __host__ __device__
#else
	#define HOST_DEVICE
#endif

#endif
