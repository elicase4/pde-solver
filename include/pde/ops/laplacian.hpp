#pragma once

#include "ops.hpp"
#include <vector>
#include <cstddef>

namespace pde{
	namespace ops{

		template<int stencil_dim>
		class Laplacian : public Operator<Laplacian<stencil_dim>>{
			public:
				Laplacian(const mesh::Mesh& mesh_): Operator<Laplacian<stencil_dim>>(mesh_) {};
				
				int stencilRadius() const{
					return stencil_dim;
				}
		};
	}
}
