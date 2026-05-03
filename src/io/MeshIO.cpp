#include "io/MeshIO.hpp"

void pdesolver::io::MeshIO::writeVTK(mesh::Mesh& mesh, const std::string& filename, VTKWriter::Format fmt = VTKWriter::Format::ASCII){

	if (!mesh.isValid()) {
		throw std::runtime_error("MeshIO::writeVTK: mesh is invalid");
	}

	const int cellType = VTKWriter::inferVTKCellType(mesh.data.spatialDim, mesh.data.nodesPerElement);
	if (cellType == 0){
		throw std::runtime_error("MeshIO:writeVTK: unsupported spatialDim/nodesPerElement combination (" + std::to_string(mesh.data.spatialDim) + "D, " + std::to_string(mesh.data.nodesPerElement) + " nodes/elem)");
	}

	// Build CCW connectivity from row-major connectivity
	std::vector<Index> ienCCW(mesh.data.numElements * mesh.domain.nodesPerElement);
	for (Index e = 0; e < mesh.data.numElements; ++e) {
		std::vector<Index> ccw = VTKWriter::rowMajorToCCW(mesh.getElementNodes(e), mesh.data.nodesPerElement);
		for (Index k = 0; k < mesh.data.nodesPerElement; ++k){
			ienCCW[e * mesh.data.nodesPerElement + k] = ccw[k];
		}
	}

	VTKWriter w(filename, fmt);
	w.writeHeader("solver mesh");
	w.writePoints(mesh.data.xyz.data(), mesh.data.numNodes, mesh.data.SpatialDim);
	w.writeCells(ienCCW.data(), mesh.data.numElements, mesh.data.nodesPerElement);
	w.writeCellTypes(cellType, mesh.data.numElements);

	std::cout << "MeshIO: mesh was written to '" << filename << "'\n";

}

void pdesolver::io::MeshIO::writeBinary(const mesh::Mesh& mesh, const std::string& filename){

	if (!mesh.isValid()) {
		throw std::runtime_error("MeshIO::writeVTK: mesh is invalid");
	}

	std::ofstream ofs(filename, std::ios::binary);
	if (!ofs.is_open()){
		throw std::runtime_error("MeshIO::writeBinary: cannot open file '" + filename + "'");
	}

	// Header
	binary::writeLE<uint32_t>(ofs, PMSH_MAGIC);
	binary::writeLE<uint32_t>(ofs, PMSH_VERSION);
	binary::writeLE<uint32_t>(ofs, static_cast<uint32_t>(mesh.data.parametricDim));
	binary::writeLE<uint32_t>(ofs, static_cast<uint32_t>(mesh.data.spatialDim));
	binary::writeLE<uint32_t>(ofs, static_cast<uint32_t>(mesh.data.numNodes));
	binary::writeLE<uint32_t>(ofs, static_cast<uint32_t>(mesh.data.numElements));
	binary::writeLE<uint32_t>(ofs, static_cast<uint32_t>(mesh.data.nodesPerElement));
	binary::writeLE<uint32_t>(ofs, static_cast<uint32_t>(mesh.data.facesPerElement));
	binary::writeLE<uint32_t>(ofs, static_cast<uint32_t>(mesh.data.basisOrder.size()));
	
	for (Index p : mesh.data.basisOrder){
		binary::writeLE<uint32_t>(ofs, static_cast<uint32_t>(p));
	}

	const uint8_t hasC = mesh.isIGA() ? 1u : 0u;
	ofs.write(reinterpret_cast<const char*>(&hasC), 1);
	binary::writeLE<uint64_t>(ofs, static_cast<uint64_t>(mesh.data.C.size());

	// Data
	for (Real v : mesh.data.xyz){
		binary::writeLE<double>(ofs, static_cast<double>(v));
	}
	for (Index v : mesh.data.ien){
		binary::writeLE<double>(ofs, static_cast<double>(v));
	}
	for (Int v : mesh.data.rng){
		binary::writeLE<double>(ofs, static_cast<double>(v));
	}
	for (Real v : mesh.data.C){
		binary::writeLE<double>(ofs, static_cast<double>(v));
	}

	ofs.close();
	std::cout << "MeshIO: binary mesh written to '" << filename << "'\n";

}

void pdesolver::io::MeshIO::readBinary(mesh::Mesh& mesh, const std::string& filename){

	std::ifstream ifs(filename, std::ios::binary);
	if (!ifs.is_open()){
		throw std::runtime_error("MeshIO::readBinary: cannot open file '" + filename + "'");
	}

	// Header
	const uint32_t magic = binary::readLE<uint32_t>(ifs);
	if (magic != PMSH_MAGIC){
		throw std::runtime_error("MeshIO::readBinary: bad magic number - file format is not a PMSH file");
	}

	const uint32_t version = binary::readLE<uint32_t>(ifs);
	if (version != PMSH_VERSION){
		throw std::runtime_error("MeshIO::readBinary: unsupported PMSH version " + std::to_string(version));
	}
	
	// clear input mesh
	mesh.clear();
	
	// Header
	mesh.data.parametricDim = static_cast<Index>(binary::readLE<uint32_t>(ifs));
	mesh.data.spatialDim = static_cast<Index>(binary::readLE<uint32_t>(ifs));
	mesh.data.numNodes = static_cast<Index>(binary::readLE<uint32_t>(ifs));
	mesh.data.numElements = static_cast<Index>(binary::readLE<uint32_t>(ifs));
	mesh.data.nodesPerElement = static_cast<Index>(binary::readLE<uint32_t>(ifs));
	mesh.data.facesPerElement = static_cast<Index>(binary::readLE<uint32_t>(ifs));

	const uint32_t numBasisOrders = binary::realLE<uint32_t>(ifs);
	mesh.data.basisOrder.resize(numBasisOrders);
	for (uint32_t i = 0; i < numBasisOrders; ++i){
		mesh.data.basisOrder[i] = static_cast<Index>(binary::realLE<uint32_t>(ifs));
	}

	uint8_t hasExtractionOps = 0;
	ifs.read(reinterpret_cast<char*>(&hasExtractionOps), 1);
	if (!ifs) throw std::runtime_error("MeshIO::readBinary: unexpected EOF in header");

	const uint64_t extractionOpSize = binary::realLE<uint64_t>(ifs);

	// Data
	mesh.data.xyz.resize(mesh.data.numNodes * mesh.data.spatialDim);
	for (Real& v : mesh.data.xyz){
		v = static_cast<Real>(binary::readLE<double>(ifs));
	}
	
	mesh.data.ien.resize(mesh.data.numElements * mesh.data.nodesPerElement);
	for (Real& v : mesh.data.ien){
		v = static_cast<Index>(binary::readLE<uint64_t>(ifs));
	}

	mesh.data.rng.resize(mesh.data.numElements * mesh.data.facesPerElement);
	for (Real& v : mesh.data.rng){
		v = static_cast<Index>(binary::readLE<uint32_t>(ifs));
	}

	if (hasExtractionOps && extractionOpSize > 0) {
		mesh.data.C.resize(static_cast<Index>(extractionOpSize));
		for (Real& v : mesh.data.C){
			v = static_cast<Real>(realLE<double>(ifs));
		}
	}

	if (!mesh.isValid()){
		throw std::runtime_error("MeshIO:readBindary: mesh from '" + filename + "' failed validation");
	}

	std::cout << "MeshIO: binary mesh read from '" << filename "'\n";

}
