#include "fem/assembly/Assembler.hpp"

namespace pdesolver::fem::assembly {

template<>
class Assembler<linalg::types::backend::CPU> {
public:

	linalg::types::CSRMatrix<Real, linalg::types::backend::CPU> createMatrix(){
		
		// allocate matrix
		linalg::types::CSRMatrix<Real, linalg::types::backend::CPU> K(topoDOF_.numFreeDOFs(), topoDOF_.numFreeDOFs());
		
		// allocate adjacency list
		std::vector<std::vector<Index>> adjList(topoDOF_.numFreeDOFs());

		// build adjacency list
		for (Index e = 0; e < mesh_.data.numElements; ++e){
			
			const Index* nodeIDs = mesh_.getElementNodes(e);
			
			for (Index i = 0; i < mesh_.data.nodesPerElement; ++i){
				for (Index j = 0; j < topoDOF_.dofsPerNode(); ++j){
					
					Index TdofIDi = topoDOF_.getNodeDOF(nodeIDs[i], j);
					Index AdofIDi;
					if (topoDOF_.isConstrained(TdofIDi)) {
						continue;
					} else {
						AdofIDi = topoDOF_.toAlgebraic(TdofIDi);
					}

					for (Index k = 0; k < mesh_.data.nodesPerElement; ++k){
						for (Index l = 0; l < topoDOF_.dofsPerNode(); ++l){
					
							Index TdofIDk = topoDOF_.getNodeDOF(nodeIDs[k], l);
							Index AdofIDk;
							if (topoDOF_.isConstrained(TdofIDk)) {
								continue;
							} else {
								AdofIDk = topoDOF_.toAlgebraic(TdofIDk);
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
		for (Index i = 0; i < topoDOF_.numFreeDOFs(); ++i){
			for (auto col : adjList[i]){
				K.colIdx()[offset++] = col;
			}
			K.rowPtr()[i+1] = offset;
		}

		return K;
	}

	linalg::types::Vector<Real, linalg::types::backend::CPU> createVector(){

		linalg::types::Vector<Real, linalg::types::backend::CPU> F(topoDOF_.numFreeDOFs());
		return F;

	}

	template<eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form, typename Quadrature>
	void assembleMatrix(const Real time, const Model& model, const Form& form, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, linalg::types::CSRMatrix<Real, linalg::types::backend::CPU>& K){
		
		// allocate local space for Ke
		linalg::types::Matrix<Real, linalg::types::backend::CPU> Ke( (EvalEle::NodesPerElement * topoDOF_.dofsPerNode()), (EvalEle::NodesPerElement * topoDOF_.dofsPerNode()) );

		// allocate local space for Ue
		linalg::types::Vector<Real, linalg::types::backend::CPU> Ue( (EvalEle::NodesPerElement * topoDOF_.dofsPerNode()) );

		// zero-out data in K
		K.zero();

		// element loop
		for (Index e = 0; e < mesh_.data.numElements; ++e){
			
			// extract node coordinates
			const Index* nodeIDs = mesh_.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * EvalEle::NodesPerElement];
			
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				
				const Real* nodeCoordsPtr = mesh_.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// zero-out Ke
			Ke.zero();
			
			// zero-out Ue
			Ue.zero();
			
			// gather U into Ue
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				for (Index j = 0; j < topoDOF_.dofsPerNode(); ++j){
					
					Index TdofIDi = topoDOF_.getNodeDOF(nodeIDs[i], j);
					if (topoDOF_.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF_.toAlgebraic(TdofIDi);
					
					Ue.data()[i*(topoDOF_.dofsPerNode()) + j] = U.data()[AdofIDi];

				}
			}

			// get element data
			EvalEle evalE;
			evalE.bindElement(nodeCoords, time);
			
			// qp data
			EvalQP qp(evalE);
			Real xi[Quadrature::NumPointsTotal*EvalEle::ParametricDim];
			Real w[Quadrature::NumPointsTotal];
			Quadrature::getPoints(xi);
			Quadrature::getWeights(w);
			
			// quadrature loop
			for (Index q = 0; q < Quadrature::NumPointsTotal; ++q){
				qp.evaluate(&xi[EvalEle::ParametricDim*q], w[q]);
				model.eval(qp);
				model.evalGradient(qp);
				form.computeElementLevel(qp, Ue.data(), Ke.data());
			}
			
			// scatter Ke into K
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				for (Index j = 0; j < topoDOF_.dofsPerNode(); ++j){
					
					Index TdofIDi = topoDOF_.getNodeDOF(nodeIDs[i], j);
					if (topoDOF_.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF_.toAlgebraic(TdofIDi);

					Index rowStart = K.rowPtr()[AdofIDi];
					Index rowEnd = K.rowPtr()[AdofIDi + 1];
					
					for (Index k = 0; k < EvalEle::NodesPerElement; ++k){
						for (Index l = 0; l < topoDOF_.dofsPerNode(); ++l){
					
							Index TdofIDk = topoDOF_.getNodeDOF(nodeIDs[k], l);
							if (topoDOF_.isConstrained(TdofIDk)) continue;
							Index AdofIDk = topoDOF_.toAlgebraic(TdofIDk);
							
							for (Index p = rowStart; p < rowEnd; ++p){
								if (K.colIdx()[p] == AdofIDk){
									K.data()[p] += Ke.data()[(i*topoDOF_.dofsPerNode() + j)*(EvalEle::NodesPerElement * topoDOF_.dofsPerNode()) + (k*topoDOF_.dofsPerNode() + l)];
									break;
								}
							}

						}
					}

				}
			}
			
		}

	}
	
	template<eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form, typename Quadrature>
	void assembleVector(const Real time, const Model& model, const Form& form, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate local space for Fe
		linalg::types::Vector<Real, linalg::types::backend::CPU> Fe( (EvalEle::NodesPerElement * topoDOF_.dofsPerNode()) );

		// allocate local space for Ue
		linalg::types::Vector<Real, linalg::types::backend::CPU> Ue( (EvalEle::NodesPerElement * topoDOF_.dofsPerNode()) );

		// zero-out data in F
		F.zero();

		// element loop
		for (Index e = 0; e < mesh_.data.numElements; ++e){
			
			// extract node coordinates
			const Index* nodeIDs = mesh_.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * EvalEle::NodesPerElement];
			
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				
				const Real* nodeCoordsPtr = mesh_.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// zero-out Ke
			Fe.zero();
			
			// zero-out Ue
			Ue.zero();
			
			// gather U into Ue
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				for (Index j = 0; j < topoDOF_.dofsPerNode(); ++j){
					
					Index TdofIDi = topoDOF_.getNodeDOF(nodeIDs[i], j);
					if (topoDOF_.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF_.toAlgebraic(TdofIDi);
					
					Ue.data()[i*(topoDOF_.dofsPerNode()) + j] = U.data()[AdofIDi];

				}
			}

			// get element data
			EvalEle evalE;
			evalE.bindElement(nodeCoords, time);
			
			// qp data
			EvalQP qp(evalE);
			Real xi[Quadrature::NumPointsTotal*EvalEle::ParametricDim];
			Real w[Quadrature::NumPointsTotal];
			Quadrature::getPoints(xi);
			Quadrature::getWeights(w);
			
			// quadrature loop
			for (Index q = 0; q < Quadrature::NumPointsTotal; ++q){
				qp.evaluate(&xi[EvalEle::ParametricDim*q], w[q]);
				model.eval(qp);
				model.evalGradient(qp);
				form.computeElementLevel(qp, Ue.data(), Fe.data());
			}
			
			// scatter Fe into F
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				for (Index j = 0; j < topoDOF_.dofsPerNode(); ++j){
					
					Index TdofIDi = topoDOF_.getNodeDOF(nodeIDs[i], j);
					if (topoDOF_.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF_.toAlgebraic(TdofIDi);
					
					F.data()[AdofIDi] += Fe.data()[i*(topoDOF_.dofsPerNode()) + j];

				}
			}
			
		}
	
	}

}; // class Assembler <linalg::types::backend::CPU>

} // namespace pdesolver::fem::assembly
