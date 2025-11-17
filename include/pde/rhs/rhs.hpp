#pragma once

#include <functional>
#include "pde/mesh/mesh.hpp"

namespace pde{
	namespace rhs {
	
		template<int DIM>
		class RHS {
			public:
				RHS(const pde::mesh::Mesh<DIM>& mesh_, const std::function<double(const std::array<double, DIM>&)> func_, double* _buffer_base): mesh(mesh_), func(func_), buffer_base(_buffer_base) {};
				
				void eval(double* f) const{
					std::ptrdiff_t flat = f - buffer_base;
					auto idx = mesh.flat2idx(static_cast<std::size_t>(flat));
					auto coord = mesh.idx2coord(idx[0], idx[1], idx[2]);
					f[0] = func(coord);
				}
			
				const pde::mesh::Mesh<DIM>& mesh;
				std::function<double(const std::array<double, DIM>&)> func;
				double* buffer_base = nullptr;
		};

	}
}
