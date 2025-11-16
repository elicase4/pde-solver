#include "pde/grid.hpp"

// cartesian to float indexing 2d
template <>
inline std::size_t pde::Grid<2>::idx(int i, int j, int) const {
	return (std::size_t(i*N[1]) + std::size_t(j));
}

// cartesian to float indexing 3d
template <>
inline std::size_t pde::Grid<3>::idx(int i, int j, int k) const {
	return (std::size_t(i*N[1]) + std::size_t(j) + std::size_t(k*N[0]*N[1]));
}

// check bounds 2d
template <>
inline bool pde::Grid<2>::in_bounds(int i, int j, int) const {
	return ((i >= 0 && i < N[0]) && (j >= 0 && j < N[1]));
}

// check bounds 3d
template <>
inline bool pde::Grid<3>::in_bounds(int i, int j, int k) const {
	return ((i >= 0 && i < N[0]) && (j >= 0 && j < N[1]) && (k >= 0 && k < N[2]));
}

// index to coordinates 2d
template <>
inline std::array<double,2> coord(int i, int j, int) const {
	std::array<double,2> out {};
	out[0] = xmin[0] + (i + 1)*h[0];
	out[1] = xmin[1] + (j + 1)*h[1];
}

// index to coordinates 3d
template <>
inline std::array<double,3> coord(int i, int j, int k) const {
	std::array<double,3> out {};
	out[0] = xmin[0] + (i + 1)*h[0];
	out[1] = xmin[1] + (j + 1)*h[1];
	out[2] = xmin[2] + (j + 1)*h[2];
}

// instantiate templates
template class pde::Grid<2>;
template class pde::Grid<3>;
