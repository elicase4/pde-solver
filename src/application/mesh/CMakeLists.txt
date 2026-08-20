cmake_minimum_required(VERSION 3.18)

find_package(yaml-cpp REQUIRED)

# Mesh application library
add_library(pde_mesh_app
    MeshConfigParser.cpp
    MeshDispatcher.cpp
    GmshMeshGenerator.cpp
)

target_include_directories(pde_mesh_app
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(pde_mesh_app
    PUBLIC
        pde_mesh
        pde_io
        pde_solver
        yaml-cpp
)

# Mesh executable
add_executable(mesh
    MeshApplication.cpp
)

target_include_directories(mesh
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(mesh
    PRIVATE
        pde_mesh_app
)

# Install target
install(TARGETS mesh
    RUNTIME DESTINATION bin
)
