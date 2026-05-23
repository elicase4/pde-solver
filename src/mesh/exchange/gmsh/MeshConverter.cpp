#include "mesh/exchange/gmsh/MeshConverter.hpp"

void pdesolver::mesh::exchange::gmsh::MeshConverter::toSolverMesh(pdesolver::mesh::Mesh& mesh, const pdesolver::mesh::exchange::gmsh::IntermediateMesh& input, const std::unordered_map<Int, Int>& physicalGroupMap) {

	if (input.empty()){
		throw std::runtime_error("MeshConverter: input IntermediateMesh is empty");
	}

	mesh.clear();
	buildConnectivity(mesh, input);
	buildBoundaryTags(mesh, input, physicalGroupMap);

	if (!mesh.isValid()){
		throw std::runtime_error("MeshConverter: produced invalid mesh");
	}

}

void pdesolver::mesh::exchange::gmsh::MeshConverter::buildConnectivity(pdesolver::mesh::Mesh& mesh, const pdesolver::mesh::exchange::gmsh::IntermediateMesh& input) {

	const ElementBlock* protoBlock = nullptr;
	
	Index totalElems = 0;
	
	// compute the total number of solver elements
	for (const auto* eb : input.cellBlocks()) {
	
		// check that cell blocks are assigned and are of uniform type
		if (!protoBlock) {
			protoBlock = eb;
		} else if (eb->type != protoBlock->type) {
			throw std::runtime_error("MeshConverter: mixed cell types are not supported");
		}

		totalElems += eb->elementIDs.size():

	}

	// set mesh metadata
	mesh.data.parametricDim = input.parametricDim;
	mesh.data.spatialDim = input.spatialDim;
	mesh.data.basisOrder = pdesolver::mesh::exchange::gmsh::basisOrder(protoBlock->type);
	mesh.data.nodesPerElement = protoBlock->nodesPerElement;
	mesh.data.facesPerElement = pdesolver::mesh::exchange::gmsh::facesPerElement(protoBlock->type);

	// set coordinates
	const Index numNodes = input.xyz.size() / 3;
	mesh.data.numNodes = numNodes;
	mesh.data.xyz.resize(numNodes * input.spatialDim);
	for (Index n = 0; n < numNodes; ++n){
		for (Index d = 0; d < input.spatialDim; ++d){
			mesh.data.xyz[n*input.spatialDim + d] = input.xyz[n*3 + d];
		}
	}

	// set connectivity
	mesh.data.numElements = totalElems;
	mesh.data.ien.resize(totalElems * protoBlock->nodesPerElement);
	Index elemOffset = 0;
	for (const auto* eb : input.cellBlocks()){

		const Index nElem = eb->elementIDs.size();
		for (Index e = 0; e < nElem; ++e){
			
			// get element connectivity
			const Index* conn = &eb->connectivity[e * eb->nodesPerElement];
			
			// reorder connectivity to row major
			auto reordered = reorderConnectivity(conn, eb->type);
			
			// store node id for each connectivity entry
			for (Index n = 0; n < eb->nodesPerElement; ++n){

				mesh.data.ien[(elemOffset + e)*mesh.data.nodesPerElement + n] = reordered[n];

			}

		}

		elemOffset += nElem;

	}
	
}

void pdesolver::mesh::exchange::gmsh::MeshConverter::buildBoundaryTags(pdesolver::mesh::Mesh& mesh, const pdesolver::mesh::exchange::gmsh::IntermediateMesh& input, const std::unordered_map<Int, Int>& physicalGroupMap) {

	const Index fpe = mesh.data.facesPerElement;
	const Index nElem = mesh.data.numElements;
	
	// Initialize rng to -1 before filling
	mesh.data.rng.assign(nElem * fpe, -1);

	// Build face-set -> (elemID, localFace) lookup from volume mesh. Key is sorted node IDs of a face
	std::unordered_map<std::vector<Index>, std::pair<Index, Index>, pdesolver::mesh::exchange::gmsh::VecHash> faceMap;
	faceMap.reserve(nElem * fpe);

	// fill in the face map
	for (Index e = 0; e < nElem; ++e){
		const Index* elemNodes = mesh.getElementNodes(e);
		for (Index f = 0; f < fpe; ++f){
			// TODO: fix this function
			auto faceNodes = pdesolver::mesh::exchange::gmsh::localFaceNodes(elemNodes, mesh.data.nodesPerElement, f);
			std::sort(faceNodes.begin(), faceNodes.end());
			faceMap[faceNodes] = {e, f};
		}
	}

	// iterate over the boundary element blocks
	for (const auto* eb : input.boundaryBlocks()) {

		// resolve solver boundary tag with input map
		Int solverTag = eb->physicalTag; // use raw physical tag by default
		
		// check input physicalGroupMap for solver boundary ID
		auto mapIt = physicalGroupMap.find(eb->physicalTag);
		if (mapIt != physicalGroupMap.end()){
			solverTag = mapIt->second;
		} else {
			solverTag = eb->physicalTag;
		}
		
		// loop over face elements in this block
		const Index nFaceElem = eb->elementIDs.size();
		for (Index fe = 0; fe < nFaceElem; ++fe){

			std::vector<Index> key(eb->nodesPerElement);
			
			// get connectivity key for each node
			for (Index n = 0; n < eb->nodesPerElement; ++n){

				key[n] = eb->connectivity[fe*eb->nodesPerElement + n];

			}
			
			// get info about face node from face map via gmsh face node id lookup
			auto it = faceMap.find(key);

			if (it == faceMap.end()){
				continue;
			}

			auto [elemID, localFace] = it->second;

			// input solver tag into rng
			mesh.data.rng[elemID*fpe + localFace] = solverTag;

		}

	}

}

std::vector<Index> pdesolver::mesh::exchange::gmsh::MeshConverter::reorderConnectivity(const Index* conn, pdesolver::mesh::exchange::gmsh::ElementType type){

	using ET = mesh::exchange::gmsh::ElementType;

	switch (type) {

		case ET::QuadP1:
			return {conn[0], conn[1], conn[3], conn[2]};


		case ET::HexP1:
			return {conn[0], conn[1], conn[3], conn[2],
					conn[4], conn[5], conn[6], conn[7]};
		
		// TODO: impelment remaining element types
		case ET::TriP1:
		case ET::TetP1:
		case ET::TriP2:
		case ET::QuadP2:
		case ET::TetP2:
		case ET::HexP2:
			return std::vector<Index>(conn, conn + nodesPerElement(type));

		default:
			return std::vector<Index>(conn, conn + nodesPerElement(type));
	
	}

}
