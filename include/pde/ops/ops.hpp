#pragma once
#include "mesh/mesh.hpp"

namespace pde {
	namespace ops {
		
		template<typename Derived>
		class Operator {
			protected:
				const mesh::Mesh& mesh;
			
			public:
				Operator(const mesh::Mesh& mesh_): mesh(mesh_) {};
				
				void apply(const double* u, double* Au) const{
					static_cast<const Derived*>(this)->apply(u, Au);
				}

				int stencilRadius() const{
					return static_cast<const Derived*>(this)->stencilRadius();
				}
		};
	}
}
