namespace pdesolver::topology {

TopologicalDOF::TopologicalDOF(const mesh::Mesh& mesh, Index dofsPerNode) : mesh_(mesh), dofsPerNode_(dofsPerNode) {
	
	numGlobalDOFs_ = mesh_.data.numNodes * dofsPerNode_;
	numFreeDOFs_ = numGlobalDOFs_;
	
	// initialize mapping
	topoToAlg_.resize(numGlobalDOFs_);
	for (Index i = 0; i < numGlobalDOFs_; ++i){
		topoToAlg_[i] = i;
	}

}

void TopologicalDOF::getElementDOFs(Index elemId, Index* dofs) const {

	const Index* nodes = mesh_.getElementNodes(elemId);
	const Index npe = mesh_.data.nodesPerElement;

	Index k = 0;
	for (Index a = 0; a < npe; ++a) {
		for (Index c = 0; c < dofsPerNode_; ++c) {
			dofs[k++] = nodes[a] * dofsPerNode_ + c;
		}
	}

}

template<typename Element>
void TopologicalDOF::buildConstraints(const fem::boundary::BoundaryRegistry& bcRegistry){

	std::set<Index> constrainedSet;
	Index faceNodes[Element::NodesPerElement];

	// loop over elements and faces to mark constrained dofs
	for (Index e = 0; mesh_.data.numElements; ++e) {
		for (Index f = 0; f < mesh_.data.facesPerElement; ++f){
			
			// check if element has boundary nodes
			if (!mesh_.isOnBoundary(e,f)) continue;

			// check for essential bcs
			Int* tagPtr = mesh_.getBoundaryTag(e);
			Int tag = tagPtr[f];
			if (!bcRegistry.isEssential(tag)) continue;

			// get face nodes
			Index npf = Element::getNodesPerFace(f);
			Element::getFaceNodes(f, faceNodes);

			// get element nodes
			const Index* elemNodes = mesh_.getElementNodes(e);

			// mark all constrained DOFs on face f
			for (Index i = 0; i < npf; ++i) {
				Index globalNode = elemNodes[faceNodes[i]];
				for (Index c = 0; c < dofsPerNode_; ++c){
					Index dof = getNodeDOF(globalNode, c);
					constrainedSet.insert(dof);
					constraintTags_[dof] = tag;
				}
			}

		}
	}

	// build algebraic numering from free dofs
	Index algIndex = 0;
	algToTopo_.clear();
	algToTopo_.reserve(numGlobalDOFs_ - constrainedSet.size());
	
	// loop over topological dofs to build algebraic mapping
	for (Index topoDOF = 0; topoDOF < numGlobalDOFs_; ++topoDOF) {
		if (constrainedSet.find(topoDOF) != constrainedSet.end()){
			topoToAlg_[topoDOF] = -1; // mark dof as constrained
		} else {
			topoToAlg_[topoDOF] = algIndex; // map to free dof
			algToTopo_.push_back(topoDOF);
			++algIndex;
		}
	}
	
	// set num of free dofs to final tally
	numFreeDOFs_ = algIndex;

}

Int TopologicalDOF::getConstraintTag(Index topoDOF) const {
	auto it = constraintTags_.find(topoDOF);
	return (it != constraintTags_.end()) ? it->second : -1;
}

} // namespace pdesolver::topology
