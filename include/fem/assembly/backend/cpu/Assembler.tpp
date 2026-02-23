namespace pdesolver::fem::assembly {

template<>
class Assembler<backend::CPU> {
public:

	template<typename EvalElement>
	static linalg::types::CSRMatrix<Real, backend::CPU> createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){
		
		// allocate matrix
		linalg::types::CSRMatrix<Real, linalg::types::backend::CPU> K(topoDOF.numFreeDOFs(), topoDOF.numFreeDOFs());
		
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

	static linalg::types::Vector<Real, backend::CPU> createOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

		linalg::types::Vector<Real, linalg::types::backend::CPU> O(topoDOF.numFreeDOFs());
		return O;

	}

	static linalg::types::Vector<Real, backend::cpu> createRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){

		linalg::types::Vector<Real, linalg::types::backend::CPU> F(topoDOF.numFreeDOFs());
		return F;

	}

	template<typename EvalElement, typename Form>
	static void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::CSRMatrix<Real, backend::CPU>& K){
		
		// allocate local space
		linalg::types::Matrix<Real, linalg::types::backend::CPU> Ke( (EvalElement::NumNodes * topoDOF.dofsPerNode()), (EvalElement::NumNodes * topoDOF.dofsPerNode()) );

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){
			
			Index* nodeIDs = mesh.getBasisNodes(e);
			Real nodeCoords[EvalElement::SpatialDim * EvalElement::NumNodes];
			
			// extract node coordinates
			for (Index i = 0; i < EvalElement::NumNodes; ++i){
				
				Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalElement::SpatialDim; ++sD){
					nodeCoords[EvalElement::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}
			
			EvalElement evalE;
			Form form;
			evalE.bindElement(nodeCoords, time);
			evalE.quadLoop<Form>(form);
			
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

	template<typename EvalElement, typename Form>
	static void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, backend::CPU>& O){

	}

	template<typename EvalElement, typename Form>
	static void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, backend::CPU>& F){

	}

}; // class Assembler <backend::CPU>

} // namespace pdesolver::fem::assembly
