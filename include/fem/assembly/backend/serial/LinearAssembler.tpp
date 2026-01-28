namespace pdesolver::fem::assembly {

template<typename Element, typename Quadrature, typename Geometry, typename EvalContext>
linalg::types::SparseMatrix<Real, Serial> LinearAssembler<Element, Quadrature, Geometry, EvalContext, Serial>::createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){
	
	// allocate matrix
	linalg::types::SparseMatrix<Real, Serial> K(topoDOF.numFreeDOFs(), topoDOF.numFreeDOFs());
	
	// allocate adjacency list
	std::vector<std::vector<Index>> adjList(topoDOF.numFreeDOFs());

	// build adjacency list
	for (Index e = 0; e < mesh.data.numElements; ++e){
		
		Index* nodeIDs = mesh.getElementNodes(e);
		
		for (Index i = 0; i < EvalContext::NumNodes; ++i){
			for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
				
				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
				Index AdofIDi;
				if (topoDOF.isConstrained(TdofIDi)) {
					continue;
				} else {
					AdofIDi = topoDOF.toAlgebraic(TdofIDi);
				}

				for (Index k = 0; k < EvalContext::NumNodes; ++k){
					for (Index l = 0; l < topoDOF.dofsPerNode(); ++l){
				
						Index TdofIDk = topoDOF.getNodeDOF(nodeIDs[k], l);
						Index AdofIDk;
						if (topoDOF.isConstrained(TdofIDk)) {
							continue;
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
	Index nnz = 0;
	for (auto& row : adjList){
		std::sort(row.begin(), row.end());
		row.erase(std::unique(row.begin(), row.end()), row.end());
		nnz += row.size();
	}

	// allocate remaining space
	K.resize(nnz);

	// fill sparse matrix
	Index offset = 0;
	K.rowPtr()[0] = 0;
	for (Index i = 0; i < topoDOF.numFreeDOFs(); ++i){
		for (auto col : adjList[i]){
			K.colIdx()[offset++] = col;
		}
		K.rowPtr()[i+1] = offset;
	}

	return K;
}

template<typename Element, typename Quadrature, typename Geometry, typename EvalContext>
linalg::types::Vector<Real, Serial> LinearAssembler<Element, Quadrature, Geometry, EvalContext, Serial>::createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

	linalg::types::Vector<Real, Serial> O(topoDOF.numFreeDOFs());
	return O;

}

template<typename Element, typename Quadrature, typename Geometry, typename EvalContext>
linalg::types::Vector<Real, Serial> LinearAssembler<Element, Quadrature, Geometry, EvalContext, Serial>::createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

	linalg::types::Vector<Real, Serial> F(topoDOF.numFreeDOFs());
	return F;

}

template<typename Element, typename Quadrature, typename Geometry, typename EvalContext, typename BilinearForm>
void LinearAssembler<Element, Quadrature, Geometry, EvalContext, Serial>::assembleMatrixSystem<BilinearForm>(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::SparseMatrix<Real, Serial>& K){

}

template<typename Element, typename Quadrature, typename Geometry, typename EvalContext, typename BilinearForm>
void LinearAssembler<Element, Quadrature, Geometry, EvalContext, Serial>::assembleOperatorSystem<BilinearForm>(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Serial>& O){

}

template<typename Element, typename Quadrature, typename Geometry, typename EvalContext, typename LinearForm>
void LinearAssembler<Element, Quadrature, Geometry, EvalContext, Serial>::assembleRHSVector<LinearForm>(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Serial>& F){

}

} // namespace pdesolver::fem::assembly
