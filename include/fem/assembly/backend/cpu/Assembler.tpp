namespace pdesolver::fem::assembly {

template<>
class Assembler<linalg::types::backend::CPU> {
public:

	PDE_HOST PDE_DEVICE linalg::types::CSRMatrix<Real, linalg::types::backend::CPU> createMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF){
		
		// allocate matrix
		linalg::types::CSRMatrix<Real, linalg::types::backend::CPU> K(topoDOF.numFreeDOFs(), topoDOF.numFreeDOFs());
		
		// allocate adjacency list
		std::vector<std::vector<Index>> adjList(topoDOF.numFreeDOFs());

		// build adjacency list
		for (Index e = 0; e < mesh.data.numElements; ++e){
			
			const Index* nodeIDs = mesh.getElementNodes(e);
			
			for (Index i = 0; i < mesh.data.nodesPerElement; ++i){
				for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
					
					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					Index AdofIDi;
					if (topoDOF.isConstrained(TdofIDi)) {
						continue;
					} else {
						AdofIDi = topoDOF.toAlgebraic(TdofIDi);
					}

					for (Index k = 0; k < mesh.data.nodesPerElement; ++k){
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

	PDE_HOST PDE_DEVICE linalg::types::Vector<Real, linalg::types::backend::CPU> createOperatorSystem(const mesh::Mesh&, const topology::TopologicalDOF& topoDOF){

		linalg::types::Vector<Real, linalg::types::backend::CPU> O(topoDOF.numFreeDOFs());
		return O;

	}

	PDE_HOST PDE_DEVICE linalg::types::Vector<Real, linalg::types::backend::CPU> createRHSVector(const mesh::Mesh&, const topology::TopologicalDOF& topoDOF){

		linalg::types::Vector<Real, linalg::types::backend::CPU> F(topoDOF.numFreeDOFs());
		return F;

	}

	template<eval::EvalElement EvalEle, eval::EvalQuadraturePoint EvalQP, typename Model, typename Form, typename Quadrature>
	PDE_HOST PDE_DEVICE void assembleMatrixSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::CSRMatrix<Real, linalg::types::backend::CPU>& K){
		
		// allocate local space
		linalg::types::Matrix<Real, linalg::types::backend::CPU> Ke( (EvalEle::NodesPerElement * topoDOF.dofsPerNode()), (EvalEle::NodesPerElement * topoDOF.dofsPerNode()) );

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){
			
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * EvalEle::NodesPerElement];
			
			// extract node coordinates
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				
				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}
			
			// get element data
			EvalEle evalE;
			evalE.bindElement(nodeCoords, time);
			
			// qp data
			EvalQP qp;
			Form form;
			Model model;
			Real xi[Quadrature::NumPointsTotal*EvalEle::ParametricDim];
			Real w[Quadrature::NumPointsTotal];
			Quadrature::getPoints(xi);
			Quadrature::getWeights(w);
			
			// quadrature loop
			model.eval(qp);
			for (Index q = 0; q < Quadrature::NumPointsTotal; ++q){
				qp.evaluate(evalE.nodeCoords, &xi[EvalEle::ParametricDim*q], w[q]);
				form(qp);
			}
			
			// scatter Ke into K
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
					
					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

					Index rowStart = K.rowPtr()[AdofIDi];
					Index rowEnd = K.rowPtr()[AdofIDi + 1];
					
					for (Index k = 0; k < EvalEle::NodesPerElement; ++k){
						for (Index l = 0; l < topoDOF.dofsPerNode(); ++l){
					
							Index TdofIDk = topoDOF.getNodeDOF(nodeIDs[k], l);
							if (topoDOF.isConstrained(TdofIDk)) continue;
							Index AdofIDk = topoDOF.toAlgebraic(TdofIDk);
							
							for (Index p = rowStart; p < rowEnd; ++p){
								if (K.colIdx()[p] == AdofIDk){
									K.data()[p] += Ke.data()[(i*topoDOF.dofsPerNode() + j)*(EvalEle::NodesPerElement * topoDOF.dofsPerNode()) + (k*topoDOF.dofsPerNode() + l)];
									break;
								}
							}

						}
					}

				}
			}
			
		}

	}

	template<eval::EvalElement EvalEle, eval::EvalQuadraturePoint EvalQP, typename Model, typename Form, typename Quadrature>
	PDE_HOST PDE_DEVICE void assembleOperatorSystem(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, linalg::types::backend::CPU>& O){

	}

	template<eval::EvalElement EvalEle, eval::EvalQuadraturePoint EvalQP, typename Model, typename Form, typename Quadrature>
	PDE_HOST PDE_DEVICE void assembleRHSVector(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const Real time, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

	}

}; // class Assembler <linalg::types::backend::CPU>

} // namespace pdesolver::fem::assembly
