namespace pdesolver::fem::assembly {

template<typename Basis, typename Quadrature, typename Geometry, typename EvalElement, typename Backend>
requires std::same_as<Backend, Serial>
linalg::types::SparseMatrix<Real, Serial> Assembler<Basis, Quadrature, Geometry, EvalElement, Serial>::createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){
	
	// allocate matrix
	linalg::types::SparseMatrix<Real, Serial> K(topoDOF.numFreeDOFs(), topoDOF.numFreeDOFs());
	
	// allocate adjacency list
	std::vector<std::vector<Index>> adjList(topoDOF.numFreeDOFs());

	// build adjacency list
	for (Index e = 0; e < mesh.data.numElements; ++e){
		
		Index* nodeIDs = mesh.getBasisNodes(e);
		
		for (Index i = 0; i < EvalElement::NumNodes; ++i){
			for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
				
				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
				Index AdofIDi;
				if (topoDOF.isConstrained(TdofIDi)) {
					continue;
				} else {
					AdofIDi = topoDOF.toAlgebraic(TdofIDi);
				}

				for (Index k = 0; k < EvalElement::NumNodes; ++k){
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

template<typename Basis, typename Quadrature, typename Geometry, typename EvalElement, typename Backend>
requires std::same_as<Backend, Serial>
linalg::types::Vector<Real, Serial> Assembler<Basis, Quadrature, Geometry, EvalElement, Serial>::createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

	linalg::types::Vector<Real, Serial> O(topoDOF.numFreeDOFs());
	return O;

}

template<typename Basis, typename Quadrature, typename Geometry, typename EvalElement, typename Backend>
requires std::same_as<Backend, Serial>
linalg::types::Vector<Real, Serial> Assembler<Basis, Quadrature, Geometry, EvalElement, Serial>::createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

	linalg::types::Vector<Real, Serial> F(topoDOF.numFreeDOFs());
	return F;

}

template<typename Basis, typename Quadrature, typename Geometry, typename EvalElement, typename Backend>
template<typename Form>
requires (std::same_as<Backend, Serial> && fem::form::BilinearForm<Form, EvalElement>)
void Assembler<Basis, Quadrature, Geometry, EvalElement, Serial>::assembleMatrixSystem<Form>(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::SparseMatrix<Real, Serial>& K){
	
	// allocate local space
	linalg::types::Matrix<Real, Serial> Ke( (EvalElement::NumNodes * topoDOF.dofsPerNode()), (EvalElement::NumNodes * topoDOF.dofsPerNode()) );

	// element loop
	for (Index e = 0; e < mesh.data.numElements; ++e){
		
		Index* nodeIDs = mesh.getBasisNodes(e);
		Real nodeCoords[EvalElement::Dimension * EvalElement::NumNodes];
		
		// extract node coordinates
		for (Index i = 0; i < EvalElement::NumNodes; ++i){
			
			Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

			for (Index d = 0; d < EvalElement::Dimension; ++d){
				nodeCoords[EvalElement::Dimension*i + d] = nodeCoordsPtr[d];
			}

		}
		

		Real xi[Quadrature::NPt];
		Real w[Quadrature::NPt];
		Quadrature::getPoints(xi);
		Quadrature::getWeights(w);

		EvalElement ctx;
		ctx.bindElement(nodeCoords, time);
		
		// quadrature loop
		for (Index q = 0; q < Quadrature::NPt; ++q){
			
			ctx.evaluate(xi[q], w[q]);
			Form::ComputeElementMatrix(ctx, Ke.data());
		
		}
		
		// scatter Ke into K
		for (Index i = 0; i < EvalElement::NumNodes; ++i){
			for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
				
				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
				if (topoDOF.isConstrained(TdofIDi)) continue;
				Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

				Index rowStart = K.rowPtr()[AdofIDi];
				Index rowEnd = K.rowPtr()[AdofIDi + 1];
				
				for (Index k = 0; k < EvalElement::NumNodes; ++k){
					for (Index l = 0; l < topoDOF.dofsPerNode(); ++l){
				
						Index TdofIDk = topoDOF.getNodeDOF(nodeIDs[k], l);
						if (topoDOF.isConstrained(TdofIDk)) continue;
						Index AdofIDk = topoDOF.toAlgebraic(TdofIDk);
						
						for (Index p = rowStart; p < rowEnd; ++p){
							if (K.colIdx()[p] == AdofIDk){
								K.data()[p] += Ke[(i*topoDOF.dofsPerNode() + j)*(EvalElement::NumNodes * topoDOF.dofsPerNode()) + (k*topoDOF.dofsPerNode() + l)];
								break;
							}
						}

			}
		}
		
	}

}

template<typename Basis, typename Quadrature, typename Geometry, typename EvalElement, typename Backend>
template<typename Form>
requires (std::same_as<Backend, Serial> && fem::form::BilinearForm<Form, EvalElement>)
void Assembler<Basis, Quadrature, Geometry, EvalElement, Serial>::assembleOperatorSystem<Form>(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Serial>& O){

}

template<typename Basis, typename Quadrature, typename Geometry, typename EvalElement, typename Backend>
template<typename Form>
requires (std::same_as<Backend, Serial> && fem::form::LinearForm<Form, EvalElement>)
void Assembler<Basis, Quadrature, Geometry, EvalElement, Serial>::assembleRHSVector<Form>(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, Serial>& F){

}

} // namespace pdesolver::fem::assembly
