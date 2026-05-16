#include "io/GmshReader.hpp"

pdesolver::io::GmshReader::VersionInfo pdesolver::io::GmshReader::readMeshFormat(std::istream& is) {

	double version;
	int fileType; // 0 = ASCII, 1 = binary
	int dataSize;
	is >> version >> fileType >> dataSize;

	pdesolver::io::GmshReader::VersionInfo vi;
	vi.version = version;
	vi.format = ((fileType == 0) ? (pdesolver::io::GmshReader::Format::ASCII) : (pdesolver::io::GmshReader::Format::Binary));

	if (vi.format == pdesolver::io::GmshReader::Format::Binary) {

		std::string rest;
		std::getline(is, rest);
		int32_t probe;
		is.read(reinterpret_cast<char*>(&probe), sizeof(probe));
		if (probe != 1){
			throw std::runtime_error("GmshReader: binary endian probe mismatch");
		}

	}

	return vi;

}

void pdesolver::io::GmshReader::read(pdesolver::mesh::exchange::gmsh::IntermediateMesh& mesh, const std::string& filename) {

	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open()){
		throw std::runtime_error("GmshReader: cannot open file: " + filename);
	}

	skipToSection(file, "$MeshFormat");
	auto vi = readMeshFormat(file);

	readMSH4(file, mesh, vi.format);

}

void pdesolver::io::GmshReader::readMSH4(std::istream& is, pdesolver::mesh::exchange::gmsh::IntermediateMesh& mesh, pdesolver::io::GmshReader::VersionInfo::Format fmt) {

	std::unordered_map<Index, Index> tagToIdx;
	std::unordered_map<Int, Int> entityPhys;

	std::string line;
	while(std::getline(is, line)) {
		if (line == "$PhysicalNames"){
			readPhysicalNames(is, mesh.physicalNames);
		} else if (line == "$Entities") {
			entityPhys = readEntities(is, fmt);
		} else if (line == "$Nodes" {
			readNodes(is, mesh, tagToIdx, fmt);
		} else if (line == "$Elements") {
			readElements(is, mesh, tagToIdx, entityPhys, fmt);
		}
	}

	deduceDimensions(mesh);

}
void pdesolver::io::GmshReader::readMSH2(std::istream&, pdesolver::mesh::exchange::gmsh::IntermediateMesh&, pdesolver::io::GmshReader::VersionInfo::Format) {

	throw std::runtime_error("GmshReader: MSH2 not implemented");

}
