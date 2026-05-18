#include "mesh/exchange/gmsh/MeshConverter.hpp"

void pdesolver::mesh::exchange::gmsh::MeshConverter::toSolverMesh(pdesolver::mesh::Mesh& mesh, const pdesolver::mesh::exchange::IntermediateMesh& input, const std::unordered_map<Int, Int>& physicalGroupMap = {}) {

	if (input.empty()){
		throw std::runtime_error("MeshConverter: input IntermediateMesh is empty");
	}

	mesh.clear();
	buildConnectivity(mesh, input);
	buildBoundaryTags(mesh, input, physicalGroupMap);

}

void pdesolver::mesh::exchange::gmsh::MeshConverter::buildConnectivity(pdesolver::mesh::Mesh& mesh, const pdesolver::mesh::exchange::IntermediateMesh& input) {

	// Check volume block element types and enforce uniformity. Determine the first volume block to extract metadata.
	const ElementBlock* protoBlock = nullptr;
	for (const auto& eb : input.elementBlocks) {
		if (io::gmsh::parametricDimension(eb.type) == input.parametricDim) {
			protoBlock = &eb;
			break;
		}
	}
	
	if (!protoBlock) {
		throw std::runtime_error("MeshConverter: no volume elements found");
	}

	// Fill mesh metadata
	mesh.data.parametricDim = input.parametricDim;
	mesh.data.spatialDim = input.SpatialDim;
	mesh.data.nodesPerElement = pdesolver::io::gmsh::facesPerElement(protoBlock->type);
	mesh.data.basisOrder = pdesolver::io::gmsh::basisOrder(protoBlock->type);

	// convert node coordinated from 3D input from intermediat mesh to 2D or 3D input solver mesh
	const Index nNodes = input.xyz.size() / 3;
	const Index sDim = input.spatialDim;
	
	mesh.data.numNodes = nNodes;
	mesh.data.xyz.resize(nNodes * sDim);
	for (Index n = 0; n < nNodes; ++n){
		for (Index d = 0; d < sDim; ++d){
			mesh.data.xyz[n*sDim + d] = input.xyz[n*3 + d];
		}
	}

	// Gmsh connectivity to ien
	mesh.data.numElements = totalElems;
	mesh.data.ien.resize(totalElems * protoBlock->nodesPerElement);

	Index elemOffset = 0;
	for (const auto& eb : input.elementBlocks) {

		if (pdesolver::io::gmsh::parametricDimension(eb.type) != input.parametricDim) continue;

		for (Index e = 0; e < eb.elementIDs.size(); ++e){
			for (Index n = 0; n < eb.nodesPerElement; ++n){
				mesh.data.ien[(elemOffset + e) * eb.nodesPerElement + n] = eb.connectivity[e*eb.nodesPerElement + n]
			}
			elemOffset += eb.elementIDs.size();
		}

	}

}

void pdesolver::mesh::exchange::gmsh::MeshConverter::buildBoundaryTags(pdesolver::mesh::Mesh& mesh, const pdesolver::mesh::exchange::IntermediateMesh& input, const std::unordered_map<Int, Int>& physicalGroupMap) {

	const Index npe = mesh.data.nodesPerElement;
	const Index fpn = mesh.data.facesPerElement;
	const Index nElem = mesh.data.numElements;
	
	// Initialize rng to -1 before filling
	mesh.data.rng.assign(nElem * fpn, -1);

	// Build face-set -> (elemID, localFace) lookup from volume mesh. Key is sorted node IDs of a face
	std::unordered_map<std::vector<Index>, std::pair<Index, Index>, pdesolver::io::gmsh::VecHash> faceMap;
	faceMap.reserve(nElem * fpn);

	// fill in the face map
	for (Index e = 0; e < nElem; ++e){
		const Index* en = mesh.getElementNodes(e);
		for (Index f = 0; f < fpn; ++f){
			// TODO: fix this utility function
			// std::vector<Index> faceNodes = pdesolver::io::gmsh::localFaceNodes(en, , f);
			std::sort(faceNodes.begin(), faceNodes.end());
			faceMap[faceNodes] = {e, f};
		}
	}

	// iterate over the boundary (lower dimensional) element blocks
	for (const auto& eb : input.elementBlocks) {

		if (pdesolver::io::gmsh::parametricDimension(eb.type) >= input.parametricDimension) continue;
		if (eb.type == pdesolver::mesh::exchange::gmsh::ElementType::Unknown) continue;
		
		// resolve solver boundary tag with input map
		Int solverTag = eb.physicalTag; // use raw physical tag by default
		if (!physicalGroupMap.empty()){
			auto it = physicalGroupMap.find(eb.physicalTag);
			if (it != physicalGroupMap.end()){
				solverTag = it->second;
			}
		}

		const Index faceNpe = eb.nodesPerElement;
		const Index nFaceElems = eb.elementIDs.size();

		for (Index fe = 0; fe < nFaceElems; ++fe) {

			std::vector<Index> key(faceNpe);
			for (Index n = 0; n < faceNpe; ++n){
				key[n] = eb.connectivity[fe * faceNpe + n];
			}
			std::sort(key.begin(), key.end());

			auto it = faceMap.find(key);
			if (it == faceMap.end()) continue;

			auto [elemID, localFace] = it->second;
			mesh.data.rng[elemID * fpn + localFace] = solverTag;

		}

	}

}
