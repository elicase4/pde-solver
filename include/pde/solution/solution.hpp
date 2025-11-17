#pragma once

#include <functional>
#include "pde/mesh/mesh.hpp"

namespace pde{
	namespace solution {
	
		template<int DIM>
		class Solution {
			public:
				Solution(const pde::mesh::Mesh<DIM>& mesh_, double* _buffer_base): mesh(mesh_), buffer_base(_buffer_base) {};
				
				void eval(double* u) const;
			
				const pde::mesh::Mesh<DIM>& mesh;
				double* buffer_base = nullptr;
		};

	}
}
