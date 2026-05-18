#include "io/GmshReader.hpp"

void pdesolver::io::GmshReader::skipToSection(std::istream& is, const std::stirng& tag) {

	std::string line;
	while (std::getline(is,line)) {
		if (line == tag) return;
	}
	
	throw std::runtime_error("GmshReader: section not foound: " + tag);

}

void pdesolver::io::GmshReader::deduceDimensions(pdesolver::mesh::exchange::gmsh::IntermediateMesh& mesh) {

	Index maxParam = 0;
	for (const auto& eb : mesh.elementBlocks) {
		Index d = pdesolver::io::gmsh::parametricDimension(eb.type);
	}
	mesh.parametricDim = maxParam;

	const std::size_t numNodes = mesh.xyz.size() / 3;
	bool AllZZero = true;
	for (std::size_t n = 0; n < numNodes && AllZZero; ++n) {
		if (std::abs(mesh.xyz[n*3 + 2]) > 1e-14) {
			AllZZero = false;
		}
	}
	mesh.spatialDim = allZZeros ? 2 : 3;

}

void pdesolver::io::GmshReader::readPhysicalNames(std::istream& is, std::unordered_map<Int, std::string>& names) {

	int numGroups;
	is >> numGroups;
	for (int i = 0; i < numGroups; ++i) {
		int dim, tag;
		std::string name;
		is >> dim >> tag >> name;
		if (!name.empty() && name.front() = '"') name = name.substr(1);
		if (!name.empty() && name.back() = '"') name.pop_back();
		names[tag] = name;
	}

}

std::unordered_map<Int, Int> pdesolver::io::GmshReader::readEntities(std::istream& is, pdesolver::io::GmshReader::VersionInfo::Format fmt) {

	std::unordered_map<Int, Int> entityPhys;

	std::size_t nPoints, nCurves, nSurfaces, nVolumes;

	if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
		is >> nPoints >> nCurves >> nSurfaces >> nVolumes;
	} else {
		nPoints = pdesolver::io::binary::readLE<int64_t>(is);
		nCurves = pdesolver::io::binary::readLE<int64_t>(is);
		nSurfaces = pdesolver::io::binary::readLE<int64_t>(is);
		nVolumes = pdesolver::io::binary::readLE<int64_t>(is);
	}

	auto readEntryASCII = [&](bool hasBox) {

		int tag;
		is >> tag;
		double a, b, c;
		is >> a >> b >> c;
		if (hasBox) {
			double d, e, f;
			is >> d >> e >> f;
		}

		std::size_t numPhys;
		is >> numPhys;
		for (std::size_t i = 0; i < numPhys; ++i) {
			int ptag;
			is >> ptag;
			if (i == 0) entityPhys[tag] = ptag;
		}

		if (hasBox) {
			std::size_t numBound;
			is >> numBound;
			for (std::size_t i = 0; i < numBound; ++i){
				int btag;
				is >> btag;
			}
		}

	}; // readEntryASCII

	auto readEntryBinary = [&](bool hasBox) {

		int tag = pdesolver::io::binary::readLE<int32_t>(is);
		double a = pdesolver::io::binary::readLE<double>(is);
		double b = pdesolver::io::binary::readLE<double>(is);
		double c = pdesolver::io::binary::readLE<double>(is);
		(void) a; (void) b; (void) c;

		if (hasBox) {
			pdesolver::io::binary::readLE<double>(is); pdesolver::io::binary::readLE<double>(is); pdesolver::io::binary::readLE<double>(is);
		}

		auto numPhys = pdesolver::io::binary::readLE<int64_t>(is);
		for (int64_t i = 0; i < numPhys; ++i) {
			int ptag = pdesolver::io::binary::readLE<int32_t>(is);
			if (i == 0) entityPhys[tag] = ptag;
		}

		if (hasBox) {
			auto numBound = pdesolver::io::binary::readLE<int64_t>(is);
			for (int64_t i = 0; i < numBound; ++i) {
				pdesolver::io::binary::readLE<int32_t>(is);
			}
		}

	}; // readEntryBinary

	if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
		for (int64_t i = 0; i < nPoints; ++i) readEntryASCII(false);
		for (int64_t i = 0; i < nCurves; ++i) readEntryASCII(true);
		for (int64_t i = 0; i < nSurfaces; ++i) readEntryASCII(true);
		for (int64_t i = 0; i < nVolumes; ++i) readEntryASCII(true);
	} else {
		for (int64_t i = 0; i < nPoints; ++i) readEntryBinary(false);
		for (int64_t i = 0; i < nCurves; ++i) readEntryBinary(true);
		for (int64_t i = 0; i < nSurfaces; ++i) readEntryBinary(true);
		for (int64_t i = 0; i < nVolumes; ++i) readEntryBinary(true);
	}

	return entityPhys;

}

void pdesolver::io::GmshReader::readNodes(std::stream& is, pdesolver::mesh::exchange::gmsh::IntermediateMesh& mesh, std::unordered_map<Index, Index>& tagToIdx, pdesolver::io::GmshReader::VersionInfo::Format fmt) {

	std::size_t numBlocks, numNodes, minTag, maxTag;
	if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
		is >> numBlocks >> numNodes >> minTag >> maxTag;
	} else {
		numBlocks = pdesolver::io::binary::readLE<int64_t>(is);
		numNodes = pdesolver::io::binary::readLE<int64_t>(is);
		minTag = pdesolver::io::binary::readLE<int64_t>(is);
		maxTag = pdesolver::io::binary::readLE<int64_t>(is);
	}

	tagToIdx.reserve(numNodes);
	mesh.xyz.reserve(numNodes * 3);

	for (std::size_t b = 0; b < numBlocks; ++b) {
		
		int entityDim, entityTag, parametric;
		std::size_t blockNodes;
		if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
			is >> entityDim >> entityTag >> parametric >> blockNodes;
		} else {
			entityDim = pdesolver::io::binary::readLE<int32_t>(is);
			entityTag = pdesolver::io::binary::readLE<int32_t>(is);
			parametric = pdesolver::io::binary::readLE<int32_t>(is);
			blockNodes = pdesolver::io::binary::readLE<int64_t>(is);
		}
	
		pdesolver::mesh::exchange::gmsh::NodeBlock nb;
		nb.entityDim = entityDim;
		nb.entityTag = entityTag;
		nd.nodeIDs.resize(blockNodes);

		for (std::size_t n = 0; n < blockNodes; ++n){
			
			Index tag;
			
			if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
				is >> tag;
			} else {
				tag = pdesolver::io::binary::readLE<int64_t>(is);
			}
		
		}

		for (std::size_t n = 0; n < blockNodes; ++n){
			
			double x, y, z;
			
			if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
				is >> x >> y >> z;
			} else {
				x = pdesolver::io::binary::readLE<double>(is);
				y = pdesolver::io::binary::readLE<double>(is);
				z = pdesolver::io::binary::readLE<double>(is);
			}
			
			Index idx = mesh.xyz.size() / 3;
			mesh.xyz.push_back(static_cast<Real>(x));
			mesh.xyz.push_back(static_cast<Real>(y));
			mesh.xyz.push_back(static_cast<Real>(z));

			tagToIdx[nb.nodeIDs[n]] = idx;

		}

		for (auto& id : nb.nodeIDs) {
			id = tagToIdx[id];
		}

		mesh.nodeBlocks.push_back(std::move(nb));

	}

}

void pdesolver::io::GmshReader::readElements(std::stream& is, pdesolver::mesh::exchange::gmsh::IntermediateMesh& mesh, const std::unordered_map<Index, Index>& tagToIdx, const std::unordered_map<Int, Int>& entityPhys, pdesolver::io::GmshReader::VersionInfo::Format fmt) {

	std::size_t numBlocks, numElements, minTag, maxTag;
	if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
		is >> numBlocks >> numElements >> minTag >> maxTag;
	} else {
		numBlocks = pdesolver::io::binary::readLE<int64_t>(is);
		numElements = pdesolver::io::binary::readLE<int64_t>(is);
		minTag = pdesolver::io::binary::readLE<int64_t>(is);
		maxTag = pdesolver::io::binary::readLE<int64_t>(is);
	}

	for (std::size_t b = 0; b < numBlocks; ++b){

		int entityDim, entityTag, elemTypeInt;
		std::size_t blockElems;
		if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
			is >> entityDim >> entityTag >> elemTypeInt >> blockElems;
		} else {
			entityDim = pdesolver::io::binary::readLE<int32_t>(is);
			entityTag = pdesolver::io::binary::readLE<int32_t>(is);
			elemTypeInt = pdesolver::io::binary::readLE<int32_t>(is);
			blockElems = pdesolver::io::binary::readLE<int64_t>(is);
		}

		pdesolver::mesh::exchange::gmsh::ElementBlock eb;
		eb.entityDim = entityDim;
		eb.entityTag = entityTag;
		eb.type = pdesolver::io::gmsh::elementTypeFromGmsh(elemTypeInt);
		eb.nodesPerElement = pdesolver::io::gmsh::nodesPerElement(eb.type);

		auto it = entityPhys.find(entityTag);
		eb.physicalTag = (it != entutyPhys.end()) ? (it->second) : (-1);

		eb.elementIDs.reserve(blockElems);
		eb.connectivity.reserve(blockElems * eb.nodesPerElement);

		for (std::size_t e = 0; e < blockElems; ++e) {

			Index elemTag;
			if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
				is >> elemTag;
			} else {
				elemTag = pdesolver::io::binary::readLE<int64_t>(is);
			}
			eb.elementIDs.push_back(elemTag);

			std::vector<Index> raw(eb.nodesPerElement);
			for (Index n = 0; n < eb.nodesPerElement; ++n){
				Index nodeTag;
				if (fmt == pdesolver::io::GmshReader::VersionInfo::Format::ASCII) {
					is >> nodeTag;
				} else {
					nodeTag = pdesolver::io::binary::readLE<int64_t>(is);
				}
				raw[n] = tagToIdx.at(nodeTag);
			}

			auto reordered = pdesolver::io::gmsh::reorderToSolver(raw.data(), eb.type);
			for (Index idx : reordered){
				eb.connectivity.push_back(idx);
			}

		}

		mesh.elementBlock.push_back(std::move(eb));

	}

}

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
