namespace pdesolver::fem::assembly {

linalg::types::SparseMatrix<Real, Serial> createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){
	
	// allocate space
	linalg::types::SparseMatrix<Real, Serial> K(topoDOF.numFreeDOFs(), topoDOF.numFreeDOFs());
	std::vector<std::vector<Index>> adjList;
	adjList.resize(topoDOF.numFreeDOFs());
	
	for (Index e = 0; e < mesh.data.numElements; ++e){
		Index* nodeIDs = mesh.getElementNodes(e);
		for (Index i = 0; i < Element::NodesPerElement; ++i){
			Index nodeID = nodeIDs[i];
		}
	}

	return K;
}

linalg::types::Vector<Real, Serial> createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

	linalg::types::Vector<Real, Serial> O(topoDOF.numFreeDOFs());
	return O;

}

linalg::types::Vector<Real, Serial> createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

	linalg::types::Vector<Real, Serial> F(topoDOF.numFreeDOFs());
	return F;

}

template<typename BilinearForm>
void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::SparseMatrix<Real, Serial>& K){

}

template<typename BilinearForm>
void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Serial>& O){

}

template<typename LinearForm>
void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Serial>& F){

}

} // namespace pdesolver::fem::assembly
