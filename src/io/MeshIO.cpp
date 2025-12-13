#include "io/MeshIO.hpp"

void pdesolver::io::MeshIO::writeVTK(const pdesolver::mesh::MeshBase& mesh, const std::string& filename){

	// open file
	std::ofstream ofs(filename);
	if (!ofs.is_open()){
		throw std::runtime_error("Cannot open file " + filename);
	}

	// get spatial dimension
	if ((mesh.data.spatialDim < 2) || (mesh.data.spatialDim > 3)){
		throw std::runtime_error("Unsupported spatial dimension: " + std::to_string(mesh.data.spatialDim));
	}

	// header
	ofs << "# vtk DataFile Version 3.0\n";
	ofs << "Mesh output\n";
	ofs << "ASCII\n";
	ofs << "DATASET UNSTRUCTURED_GRID\n";
	
	// points
	ofs << "POINTS " << mesh.data.numNodes << " double\n";
	for (Index n = 0; n < mesh.data.numNodes; ++n){
		const double* xyz = (double*) mesh.getNodeCoord(n);
		ofs << xyz[0] << " " << xyz[1] << " ";
		if (mesh.data.spatialDim == 2) ofs << "0.0";
		else ofs << xyz[2];
		ofs << "\n";
	}

	// cells
	ofs << "CELLS " << mesh.data.numElements << " "
		<< mesh.data.numElements * (mesh.data.nodesPerElement + 1) << "\n";
	for (Index e = 0; e < mesh.data.numElements; ++e){
		std::vector<Index> node_vec = rowMajorToCCW(mesh.getElementNodes(e), mesh.data.nodesPerElement);
		Index* nodes = node_vec.data();
		ofs << mesh.data.nodesPerElement;
		for (Index i = 0; i < mesh.data.nodesPerElement; ++i) {
			ofs << " " << nodes[i];
		}
		ofs << "\n";
	}

	// cell types
	ofs << "CELL_TYPES " << mesh.data.numElements << "\n";
	for (Index e = 0; e < mesh.data.numElements; ++e) {
		if (mesh.data.spatialDim == 2){
			switch (mesh.data.nodesPerElement){
				case 4:
					ofs << "9\n"; // quad
					break;
				default:
					ofs << "0\n"; // not found
					break;
			}
		} else if (mesh.data.spatialDim == 3){
			switch (mesh.data.nodesPerElement){
				case 8:
					ofs << "12\n"; // hex
					break;
				default:
					ofs << "0\n"; // not found
					break;
			}

		}
	}
	
	// close file
	ofs.close();
	std::cout << "Mesh written to " << filename << "\n";
}

std::vector<Index> pdesolver::io::MeshIO::rowMajorToCCW(const Index* row_major_ordering, Index num_nodes){
	switch (num_nodes){
		case (4):
			return { row_major_ordering[0], row_major_ordering[1], row_major_ordering[3], row_major_ordering[2] };
		case (8):
			return { row_major_ordering[0], row_major_ordering[1], row_major_ordering[3], row_major_ordering[2], row_major_ordering[4], row_major_ordering[5], row_major_ordering[7], row_major_ordering[6] };
		default:
			throw std::runtime_error("Unsupported element type.");
	}
}
