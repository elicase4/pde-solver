#include "pde/diffop/laplacian.hpp"

template <>
void pde::diffop::Laplacian<1,2>::apply_inst(const double* u, double* Au) const{
	Au[mesh.idx2flat(0,0)] = (1/(mesh.h[0]*mesh.h[0])) * (u[mesh.idx2flat(1,0)] - 2*u[mesh.idx2flat(0,0)] + u[mesh.idx2flat(-1,0)])
						+ (1/(mesh.h[1]*mesh.h[1])) * (u[mesh.idx2flat(0,1)] - 2*u[mesh.idx2flat(0,0)] + u[mesh.idx2flat(0,-1)]);
}

template <>
void pde::diffop::Laplacian<1,3>::apply_inst(const double* u, double* Au) const{
	Au[mesh.idx2flat(0,0,0)] = (1/(mesh.h[0]*mesh.h[0])) * (u[mesh.idx2flat(1,0,0)] - 2*u[mesh.idx2flat(0,0,0)] + u[mesh.idx2flat(-1,0,0)])
						+ (1/(mesh.h[1]*mesh.h[1])) * (u[mesh.idx2flat(0,1,0)] - 2*u[mesh.idx2flat(0,0,0)] + u[mesh.idx2flat(0,-1,0)])
						+ (1/(mesh.h[2]*mesh.h[2])) * (u[mesh.idx2flat(0,0,1)] - 2*u[mesh.idx2flat(0,0,0)] + u[mesh.idx2flat(0,0,-1)]);
}

template class pde::diffop::Laplacian<1,2>;
template class pde::diffop::Laplacian<1,3>;
