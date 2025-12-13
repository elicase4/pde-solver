#include "core/Types.hpp"
#include "mesh/BlockMesh2D.hpp"
#include "io/MeshIO.hpp"

using namespace pdesolver;

int main(int argc, char** argv){
	
	// get filename from cli
	std::string filename;
	if (argc == 2){
		filename = argv[1];
	} else {
		std::cout << "Default filename: meshblock2d.vtk will be used.\n";
		filename = "meshblock2d.vtk";
	}
	
	// initialize mesh parameters, add cli with boost later
	Index nx = 2; Index ny = 2;
	Real x0 = -1.0; Real x1 = 1.0; Real y0 = -1.0; Real y1 = 1.0;
	Index px = 1; Index py = 1;
	
	// construct blockmesh2d object
	auto mesh2d = mesh::BlockMesh2D(nx, ny, x0, x1, y0, y1, px, py);
	
	// build blockmesh2d
	mesh2d.initializeData();
	mesh2d.generateNodes();
	mesh2d.generateElements();
	mesh2d.generateBoundaryTags();
	
	std::cout << mesh2d.data.nodesPerElement << "\n";

	// write file in vtk format
	auto meshio = io::MeshIO();
	meshio.writeVTK(mesh2d, filename);

	return 0;
}
