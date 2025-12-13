add_library(core
    pde/mesh/mesh.cpp
    pde/diffop/laplacian.cpp
)

target_include_directories(core PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(solvers
    solver/jacobi.cpp
    solver/gauss_seidel.cpp
    solver/conjugate_gradient.cpp
    solver/sor.cpp
)

target_link_libraries(solvers PUBLIC core)

if (ENABLE_OMP)
    add_library(openmp_backend
        backend/openmp/apply_operator.cpp
        backend/openmp/compute_residual.cpp
        backend/openmp/update_solution.cpp
    )
    target_link_libraries(openmp_backend PUBLIC solvers)
endif()

if (ENABLE_MPI)
    add_library(mpi_backend
        backend/mpi/apply_operator.cpp
        backend/mpi/compute_residual.cpp
        backend/mpi/update_solution.cpp
        backend/mpi/decompose.cpp
        backend/mpi/halo_exchange.cpp
    )
    target_link_libraries(mpi_backend PUBLIC solvers)
endif()

if (ENABLE_CUDA)
    add_library(cuda_backend
        backend/cuda/apply_operator.cpp
        backend/cuda/compute_residual.cpp
        backend/cuda/update_solution.cpp
    )
    target_link_libraries(cuda_backend PUBLIC solvers)
endif()
