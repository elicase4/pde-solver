#pragma once

#include "pde/diffop/diffop.hpp"

namespace pde{
	namespace diffop{

		template<int stencil_rad, int DIM>
		class Laplacian : public DiffOp<Laplacian<stencil_rad, DIM>>{
			public:
				Laplacian(const mesh::Mesh<DIM>& mesh_): mesh(mesh_) {};
				
				void apply_inst(const double* u, double* Au) const;
				
				int stencilRadius_inst() const{
					return stencil_rad;
				}
			private:
				const mesh::Mesh<DIM>& mesh;
		};

	}
}
