#include "pde/rhs/rhs.hpp"

template<>
void pde::rhs::RHS<2>::eval(double* f) const{
	std::ptrdiff_t flat = f - buffer_base;
	std::array<int,2> idx = mesh.flat2idx(static_cast<std::size_t>(flat));
	std::array<double,2> coord = mesh.idx2coord(idx[0], idx[1]);
	f[0] = func(coord);
}

template<>
void pde::rhs::RHS<3>::eval(double* f) const{
	std::ptrdiff_t flat = f - buffer_base;
	std::array<int,3> idx = mesh.flat2idx(static_cast<std::size_t>(flat));
	std::array<double,3> coord = mesh.idx2coord(idx[0], idx[1], idx[2]);
	f[0] = func(coord);
}

template class pde::rhs::RHS<2>;
template class pde::rhs::RHS<3>;
