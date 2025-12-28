add_executable(poisson_serial
	poisson/poisson.exe
)

target_link_libraries(poisson_serial
	PRIVATE
		pde_mesh
		pde_io
		pde_fem
		pde_linalg
)
