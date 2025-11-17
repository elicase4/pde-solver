#include "pde/mesh/mesh.hpp"

template <>
std::size_t pde::mesh::Mesh<2>::idx2flat(int i, int j, int) const {
	return static_cast<std::size_t>(i*N[1] + j);
}

template <>
std::size_t pde::mesh::Mesh<3>::idx2flat(int i, int j, int k) const {
	return static_cast<std::size_t>(i*N[1] + j + k*N[0]*N[1]);
}

template <>
std::array<int,2> pde::mesh::Mesh<2>::flat2idx(std::size_t flat_idx) const {
	std::array<int,2> out {};
	out[1] = static_cast<int>(flat_idx % N[1]);
	out[0] = static_cast<int>((flat_idx - out[1]) / N[1]);
	return out;
}

template <>
std::array<int,3> pde::mesh::Mesh<3>::flat2idx(std::size_t flat_idx) const {
	std::array<int,3> out {};
	out[1] = static_cast<int>(flat_idx % N[1]);
	out[0] = static_cast<int>((flat_idx % N[0] - out[1]) / N[1]);
	out[2] = static_cast<int>((flat_idx - out[0] * N[1] - out[1]) / (N[1]*N[0]));
	return out;
}

template <>
bool pde::mesh::Mesh<2>::in_bounds(int i, int j, int) const {
	return ((i >= 0 && i < N[0]) && (j >= 0 && j < N[1]));
}

template <>
bool pde::mesh::Mesh<3>::in_bounds(int i, int j, int k) const {
	return ((i >= 0 && i < N[0]) && (j >= 0 && j < N[1]) && (k >= 0 && k < N[2]));
}

template <>
std::array<double,2> pde::mesh::Mesh<2>::idx2coord(int i, int j, int) const {
	std::array<double,2> out {};
	out[0] = xmin[0] + (i + 1)*h[0];
	out[1] = xmin[1] + (j + 1)*h[1];
	return out;
}

template <>
std::array<double,3> pde::mesh::Mesh<3>::idx2coord(int i, int j, int k) const {
	std::array<double,3> out {};
	out[0] = xmin[0] + (i + 1)*h[0];
	out[1] = xmin[1] + (j + 1)*h[1];
	out[2] = xmin[2] + (k + 1)*h[2];
	return out;
}

template <>
std::size_t pde::mesh::Mesh<2>::size() const {
	return static_cast<std::size_t>(N[0] * N[1]);
}

template <>
std::size_t pde::mesh::Mesh<3>::size() const {
	return static_cast<std::size_t>(N[0] * N[1] * N[2]);
}

template class pde::mesh::Mesh<2>;
template class pde::mesh::Mesh<3>;
