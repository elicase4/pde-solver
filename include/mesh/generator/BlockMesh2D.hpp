#ifndef PDESOLVER_BLOCKMESH2D_HPP
#define PDESOLVER_BLOCKMESH2D_HPP

#include "mesh/Mesh.hpp"

namespace pdesolver {
	
	namespace mesh {
		
		namespace generator {

			class BlockMesh2D : public Mesh {
			public:
				BlockMesh2D(Index nx_, Index ny_, Real x0_, Real x1_, Real y0_, Real y1_, Index px_, Index py_): nx(nx_), ny(ny_), x0(x0_), x1(x1_), y0(y0_), y1(y1_), px(px_), py(py_) {};
				
				void initializeData();

				void generateNodes();

				void generateElements();

				void generateBoundaryTags();
			
			private: 
				
				Index nx;
				Index ny;
				
				Real x0;
				Real x1;
				Real y0;
				Real y1;
				
				Index px;
				Index py;

			}; // class BlockMesh2D

		} // namespace generator

	} // namespace mesh

} // namespace pdesolver

#endif
