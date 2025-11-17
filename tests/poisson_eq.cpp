#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "pde/mesh/mesh.hpp"
#include "pde/solution/solution.hpp"
#include "pde/rhs/rhs.hpp"
#include "pde/diffop/laplacian.hpp"
#include "pde/boundary/dirichlet.hpp"

void check2DCase(int N){
	
	using Mesh2 = pde::mesh::Mesh<2>;

    Mesh2 mesh(
        {N, N},
        {-1.0,-1.0},
        {1.0,1.0}
    );

    const std::size_t n = mesh.size();
    std::vector<double> u(n);
    std::vector<double> f(n);
    std::vector<double> f_comp(n);
    std::vector<double> Lu(n);
	
	// test rhs function
	auto rhs_f = [](const std::array<double,2>& x){
		return std::sin(M_PI*x[0]) * std::sin(M_PI*x[1]);
	};

	// test dirichlet bc
	auto const_dbc_f = [](const std::array<double,2>& x){
		return 10.0;
	};
	
	// construct pde objects
	pde::diffop::Laplacian<1,2> lap(mesh);
    pde::rhs::RHS<2> f_rhs(mesh, rhs_f, &f[0]);
    pde::solution::Solution<2> u_sol(mesh, &u[0]);
	pde::boundary::Dirichlet<2> dbc(mesh, u_sol, f_rhs);
	
	// set bcs
	dbc.set(pde::mesh::BoundaryFace::XMin, const_dbc_f);
	dbc.set(pde::mesh::BoundaryFace::XMax, const_dbc_f);
	dbc.set(pde::mesh::BoundaryFace::YMin, const_dbc_f);
	dbc.set(pde::mesh::BoundaryFace::YMax, const_dbc_f);

    // test function u(x,y) = sin(pi*x) sin(pi*y)
    for (int i=0; i<mesh.N[0]; i++){
		for (int j=0; j<mesh.N[1]; j++){
			auto c = mesh.idx2coord(i,j);
			double x = c[0], y = c[1];
			
			auto face = mesh.get_boundary_face(i,j);
			if(face != pde::mesh::BoundaryFace::None){
				dbc.eval(&u[mesh.idx2flat(i,j)], &f[mesh.idx2flat(i,j)]);
			} else {
				u[mesh.idx2flat(i,j)] = (-1 / (2 * M_PI * M_PI)) * std::sin(M_PI*x)*std::sin(M_PI*y) + 10.0;
			}
			
			f_comp[mesh.idx2flat(i,j)] = std::sin(M_PI*x)*std::sin(M_PI*y);
		}
    }

    for (int i=1; i<mesh.N[0]-1; i++){
		for (int j=1; j<mesh.N[1]-1; j++){
			lap.apply(&u[mesh.idx2flat(i,j)], &Lu[mesh.idx2flat(i,j)]);
			f_rhs.eval(&f[mesh.idx2flat(i,j)]);
		}
    }

    std::cout << "grid output:\n";
    std::cout << std::setw(8) << "i"
              << std::setw(8) << "j"
              << std::setw(20) << "Numerical Laplacian"
              << std::setw(20) << "Exact Laplacian"
              << std::setw(20) << "RHS"
              << std::setw(20) << "Solution Vector"
              << "\n";
    for (int i=0; i<mesh.N[0]-1; i+=10) {
		for (int j=0; j<mesh.N[0]-1; j+=10) {
			std::size_t id = mesh.idx2flat(i,j);

			// Print values aligned
			std::cout << std::setw(8) << i
					  << std::setw(8) << j
					  << std::setw(20) << std::fixed << std::setprecision(6) << Lu[id]
					  << std::setw(20) << std::fixed << std::setprecision(6) << f_comp[id]
					  << std::setw(20) << std::fixed << std::setprecision(6) << f[id]
					  << std::setw(20) << std::fixed << std::setprecision(6) << u[id]
					  << "\n";
		}
    }
}

int main() {
	const int N = 100;
	check2DCase(N);
	
	return 0;
}    
