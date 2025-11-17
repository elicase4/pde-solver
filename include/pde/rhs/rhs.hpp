#pragma once

#include <functional>
#include "pde/mesh/mesh.hpp"

namespace pde{
	namespace rhs {
	
		template<int DIM>
		class RHS {
			public:
				RHS(const pde::mesh::Mesh<DIM>& mesh_, const std::function<double(const std::array<double, DIM>&)> func_, double* _buffer_base): mesh(mesh_), func(func_), buffer_base(_buffer_base) {};
				
				void eval(double* f) const;
			
				pde::mesh::Mesh<DIM> mesh;
				std::function<double(const std::array<double, DIM>&)> func;
				double* buffer_base;
		};

	}
}
