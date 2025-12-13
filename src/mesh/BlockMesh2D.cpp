#include "mesh/BlockMesh2D.hpp"

void pdesolver::mesh::BlockMesh2D::initializeData(){
	
	// dimensions
	data.parametricDim = 2;
	data.spatialDim = 2;

	// basis function info
	data.basisOrder = {px, py};
	data.nodesPerElement = (px + 1)*(py + 1);

	// mesh size info
	data.numNodes = (nx + 1) * (ny + 1);
	data.numElements = nx * ny;
	
	// boundary info
	data.facesPerElement = 4;

	// size array
	data.xyz.resize(data.numNodes * data.spatialDim);
	data.ien.resize(data.numElements * data.nodesPerElement);
	data.rng.resize(data.numElements * data.facesPerElement);

}

void pdesolver::mesh::BlockMesh2D::generateNodes(){
	
	// element spacing
	Real dx = (x1 - x0) / ( (Real) nx);
	Real dy = (y1 - y0) / ( (Real) ny);

	// populate xyz array
	Index node = 0;
	
	for (Index j = 0; j <= ny; ++j) {
		for (Index i = 0; i <= nx; ++i) {
			
			data.xyz[data.spatialDim * node] = x0 + ( (Real) i) * dx;
			data.xyz[data.spatialDim * node + 1] = y0 + ( (Real) j) * dy;
			
			node++;
		}
	}

}

void pdesolver::mesh::BlockMesh2D::generateElements(){
	
	// populate ien array
	for (Index ele_y = 0; ele_y < ny; ++ele_y){
		for (Index ele_x = 0; ele_x < nx; ++ele_x){
			
			Index ele = ele_y * nx + ele_x;
			
			for (Index a = 0; a < data.nodesPerElement; ++a){
				
				Index a_x  = a % (px + 1);
				Index a_y  = a / (px + 1);

				Index node_x = ele_x * px + a_x;
				Index node_y = ele_y * py + a_y;
				
				data.ien[ele * data.nodesPerElement + a] = node_y * (nx * px + 1) + node_x;
			
			}
		
		}
	}

}

void pdesolver::mesh::BlockMesh2D::generateBoundaryTags(){
	
	// boundary tags
	static constexpr Int LEFT   = 0;
	static constexpr Int RIGHT  = 1;
	static constexpr Int BOTTOM = 2;
	static constexpr Int TOP    = 3;

	// populate rng array
	for (Index ele_y = 0; ele_y < ny; ++ele_y){
		for (Index ele_x = 0; ele_x < nx; ++ele_x){
			
			// get flat id
			Index ele = ele_y * nx + ele_x;
			
			// set boundary tags
			data.rng[ele * data.facesPerElement + LEFT]   = (ele_x == 0)        ? LEFT   : -1;
			data.rng[ele * data.facesPerElement + RIGHT]  = (ele_x == (nx - 1)) ? RIGHT  : -1;
			data.rng[ele * data.facesPerElement + BOTTOM] = (ele_y == 0)        ? BOTTOM : -1;
			data.rng[ele * data.facesPerElement + TOP]    = (ele_y == (ny - 1)) ? TOP    : -1;

		}
	}
}
