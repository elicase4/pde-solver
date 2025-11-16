#include "pde/ops/laplacian.hpp"

void pde::ops::Laplacian<5>::apply(const double* u, double* Au) const{
	Au[pde::mesh::idx(0,0)] = (1/(h[0]*h[0])) * (u[pde::mesh::idx(1,0)] - 2*u[pde::mesh::idx(0,0)] + u[pde::mesh::idx(-1,0)])
							+ (1/(h[1]*h[1])) * (u[pde::mesh::idx(0,1)] - 2*u[pde::mesh::idx(0,0)] + u[pde::mesh::idx(0,-1)]);
}
