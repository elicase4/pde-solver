#include "io/FieldIO.hpp"

void pdesolver::io::FieldIO::writeVTK(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const fem::boundary::BoundaryRegistry& bcRegistry, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::string& filename, VTKWriter::Format fmt = VTKWriter::Format::ASCII) {

	if (!mesh.isValid()){
		throw std::runtime_error("FieldIO::writeVTK: mesh is invalid");
	}
	
	if (dofNames.size() != topoDOF.dofsPerNode()){
		throw std::runtime_error("FieldIO::writeVTK: dofNames.size() (" + std::to_string(dofNames.size()) + ") must equal dofsPerNode (" + std::to_string(topoDOF.dofsPerNode()) + ")");
	}

	const int cellType = VTKWriter::inferVTKCellType(mesh.data.spatialDim, mesh.data.nodesPerElement);
	if (cellType == 0){
		throw std::runtime_error("FieldIO::writeVTK: unsupported spatialDim/nodesPerElement "" combination (" + std::to_string(mesh.data.spatialDim) + "D, " + std::to_string(mesh.data.nodesPerElement) + " nodes/elem");
	}

	// get ccw connectivity
	std::vector<Index> ienCCW(mesh.data.numElements * mesh.data.nodesPerElement);
	for (Index e = 0; e < mesh.data.numElements; ++e){
		
		std::vector<Index> ccw = VTKWriter::rowMajorToCCW(mesh.getElementNodes(e), mesh.data.nodesPerElement);

		for (Index k = 0; k < mesh.data.nodesPerElement; ++k){
			ienCCW[e * mesh.data.nodesPerElement + k] = ccw[k]
		}

	}

	// reconstruct full nodal field
	const std::vector<Real> nodalField = reconstructNodalField(mesh, topoDOF, bcRegistry, time, algField);

	// write file
	VTKWriter w(filename, fmt);
	w.writeHeader("solver field output");
	w.writePoints(mesh.data.xyz.data(), mesh.data.numNodes, mesh.data.spatialDim);
	w.writeCells(ienCCW.data(), mesh.data.numElements, mesh.data.nodesPerElement);
	w.writeCellTypes(cellType, mesh.data.numElements);

	// write point data
	w.beginPointData(mesh.data.numNodes);

	// extract interleaved data from nodalField for stride-1 buffer to input to writeScalar
	std::vector<Real> componentBuf(mesh.data.numNodes);
	for (Index c = 0; c < topoDOF.dofsPerNode(); ++c) {
		for (Index n = 0; n < mesh.data.numNodes; ++n){
			componentBuf[n] = nodalField[n*topoDOF.dofsPerNode() + c];
		}
		w.writeScalar(dofNames[c], componentBuf.data(), mesh.data.numNodes);
	}

	w.endPointData();

	const char* ordStr = (topoDOF.ordering() == fem::dof::DOFOrdering::Interleaved) ? "Interleaved" : "Block";

	std::cout << "FieldIO: '" << filename << "' written - " << topoDOF.dofsPerNode() << " field(s), " << mesh.data.numNodes << " nodes, " << ordStr << " ordering\n";

}

std::vector<Real> pdesolver::io::FieldIO::reconstructNodalField(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const fem::boundary::BoundaryRegistry& bcRegistry, Real time, const Real* algField) const {

	// initialize containers for eval
	std::vector<Real> nodalField(mesh.data.numNodes * topoDOF.dofsPerNode(), Real(0.0));
	std::vector<Real> bcVal(topoDOF.dofsPerNode(), Real(0.0));

	// write free DOFs
	for (Index i = 0; i < topoDOF.numFreeDOFs(); ++i) {
		
		const Index topoDOFIdx = topoDOF.toTopological(i);
		const Index nodeId = topoDOF.getDOFnode(topoDOFIdx);
		const Index component = topoDOFIdx - (nodeId * topoDOF.dofsPerNode());

		nodalField[nodeID * topoDOF.dofsPerNode() + component] = algField[i];
	
	}

	// write constrained DOFs
	for (Index topoDOFIdx = 0; topoDOFIdx < topoDOF.numGlobalDOFs(); ++topoDOFIdx) {

		if (!topoDOF.isConstrained(topoDOFIdx)) continue;

		const Index nodeId = topoDOF.getDOFNode(topoDOFIdx);
		const Index component = topoDOFIdx - (nodeId * topoDOF.dofsPerNode());
		const Int tag = topoDOF.getconstraintTag(topoDOFIdx);
		const Real* xyz = mesh.getNodeCoord(nodeId);

		for (const auto& entry : bcRegistry.entries()) {
			
			if (entry->tag() != tag) continue;
			if (entry->componentType(component) != fem::boundary::BCCategory::Essential) continue;

			entry->eval(time, xyz, bcVal.data());
			nodalField[nodeID * topoDOF.dofsPerNode() + component] = bcVal[component];

			break;

		}

	}

	return nodalField;

}
