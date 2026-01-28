namespace pdesolver::fem::assembly {

linalg::types::SparseMatrix<Real, Serial> createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){
	
	// allocate matrix
	linalg::types::SparseMatrix<Real, Serial> K(topoDOF.numFreeDOFs(), topoDOF.numFreeDOFs());
	
	// allocate adjacency list
	std::vector<std::vector<Index>> adjList;
	adjList.resize(topoDOF.numFreeDOFs());
	for (Index i = 0; i < topoDOF.numFreeDOFs(); ++i){
		adjList[i].resize(EvalContext::NumNodes);
	}

	// compute nnz
	Index nnz;
	
	// build adjacency list
	for (Index e = 0; e < mesh.data.numElements; ++e){
		
		Index* nodeIDs = mesh.getElementNodes(e);
		
		for (Index i = 0; i < EvalContext::NumNodes; ++i){
			for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
				
				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
				Index AdofIDi;
				if (topoDOF.isConstrained(TdofIDi)) {
					break;
				} else {
					AdofIDi = toAlgebraic(TdofIDi);
				}

				for (Index k = 0; k < EvalContext::NumNodes; ++k){
					for (Index l = 0; l < topoDOF.dofsPerNode(); ++l){
				
						Index TdofIDk = topoDOF.getNodeDOF(nodeIDs[k], l);
						Index AdofIDk;
						if (topoDOF.isConstrained(TdofIDk)) {
							break;
						} else {
							AdofIDk = toAlgebraic(TdofIDk);
							adjList[AdofIDi].push_back(AdofIDk);
						}

					}
				}

			}
		}
	
	}

	// reduce adjacency list to unique pairs
	for (auto& row : adjList){
		std::sort(row.begin(), row.end());
		row.erase(std::unique(row.begin(), row.end()), row.end());
		nnz += row.size();
	}

	// allocate remaining space
	K.resize(nnz);

	// fill sparse matrix
	K.rowPtr[0] = 0;
	Index count = 0;
	for (Index i = 0; i < topoDOF.numFreeDOFs(); ++i){
		K.rowPtr[i+1] = K.rowPtr[i] + adjList[i].size();
		for (auto val : adjList[i]){
			K.colIdx[count++] = val;
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
