namespace pdesolver::fem::assembly {

template<typename Basis, typename Quadrature, typename Geometry, typename EvalContext>
linalg::types::SparseMatrix<Real, Serial> LinearAssembler<Basis, Quadrature, Geometry, EvalContext, Serial>::createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){
	
	// allocate matrix
	linalg::types::SparseMatrix<Real, Serial> K(topoDOF.numFreeDOFs(), topoDOF.numFreeDOFs());
	
	// allocate adjacency list
	std::vector<std::vector<Index>> adjList(topoDOF.numFreeDOFs());

	// build adjacency list
	for (Index e = 0; e < mesh.data.numElements; ++e){
		
		Index* nodeIDs = mesh.getBasisNodes(e);
		
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
							AdofIDk = topoDOF.toAlgebraic(TdofIDk);
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

template<typename Basis, typename Quadrature, typename Geometry, typename EvalContext>
linalg::types::Vector<Real, Serial> LinearAssembler<Basis, Quadrature, Geometry, EvalContext, Serial>::createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

	linalg::types::Vector<Real, Serial> O(topoDOF.numFreeDOFs());
	return O;

}

template<typename Basis, typename Quadrature, typename Geometry, typename EvalContext>
linalg::types::Vector<Real, Serial> LinearAssembler<Basis, Quadrature, Geometry, EvalContext, Serial>::createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

	linalg::types::Vector<Real, Serial> F(topoDOF.numFreeDOFs());
	return F;

}

template<typename Basis, typename Quadrature, typename Geometry, typename EvalContext>
template<typename Form>
void LinearAssembler<Basis, Quadrature, Geometry, EvalContext, Serial>::assembleMatrixSystem<Form>(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::SparseMatrix<Real, Serial>& K){
	
	// allocate local space
	linalg::types::Matrix<Real, Serial> Ke((EvalContext::NumNodes * topoDOF.dofsPerNode() ) * (EvalContext::NumNodes * topoDOF.dofsPerNode() ) );

	// element loop
	for (Index e = 0; e < mesh.data.numElements; ++e){
		
		Index* nodeIDs = mesh.getBasisNodes(e);
		Real* nodeCoords;
		
		// extract node coordinates
		for (Index i = 0; i < EvalContext::NumNodes; ++i){
			
			Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

			for (Index d = 0; d < EvalContext::Dimension; ++d){
				nodeCoords[EvalContext::Dimension*i + d] = nodeCoordsPtr[d];
			}

		}
		
		EvalContext::bindElement(nodeCoords, time);

		Real* xi;
		Real* w;
		Quadrature::getPoints(xi);
		Quadrature::getWeights(w);
		
		// quadrature loop
		for (Index q = 0; q < Quadrature::NPt; ++q){
			
			EvalContext::evaluate(xi[q], w[q]);
			Form::ComputeElementMatrix(EvalContext, Ke.data());
		
		}

		// scatter Ke into K
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
							AdofIDk = topoDOF.toAlgebraic(TdofIDk);
							Index row_idx = K.rowPtr()[AdofIDi];
							Index flat_idx = K.colIdx()[row_idx + AdofIDj];
							K.data()[flat_idx] = Ke[(i*topoDOF.dofsPerNode() + j)*(EvalContext::NumNodes * topoDOF.dofsPerNode()) + (k*topoDOF.dofsPerNode() + l)];
						}

			}
		}
		
	}

}

template<typename Basis, typename Quadrature, typename Geometry, typename EvalContext>
template<typename Form>
void LinearAssembler<Basis, Quadrature, Geometry, EvalContext, Serial>::assembleOperatorSystem<Form>(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const real time, linalg::types::Vector<Real, Serial>& O){

}

template<typename Basis, typename Quadrature, typename Geometry, typename EvalContext>
template<typename Form>
void LinearAssembler<Basis, Quadrature, Geometry, EvalContext, Serial>::assembleRHSVector<Form>(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const real time, linalg::types::Vector<Real, Serial>& F){

}

} // namespace pdesolver::fem::assembly
