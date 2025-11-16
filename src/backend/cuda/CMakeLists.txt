enable_language(CUDA)

set(CUDA_SOURCES
	apply_operator.cu
	update_solution.cu
	compute_residual.cu
)

add_library(cuda_backend ${CUDA_SOURCES})
target_include_directories(cuda_backend PUBLIC ${PROJECT_SOURCE_DIR}/include)
set_target_properties(cuda_backend PROPERTIES
	CUDA_SEPARABLE_OCMPILATION ON
	CUDA_STANDARD 17
)
target_link_libraries(cuda_backend PUBLIC core)
