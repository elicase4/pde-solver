namespace pdesolver::fem::boundary {

template<>
class BoundaryApplicator<linalg::types::backend::CPU>
public:
	
	template<eval::EvalElement EvalEle, typename EvalQP, typename Form, typename Basis, typename Model>
	void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& bcRegistry, const Real time, const Model& model, const Form&form, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

	}

	template<eval::EvalElement EvalEle, typename EvalQP, typename Form, typename Basis, typename Quadrature>
	void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& bcRegistry, const Real time, const Form& form, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate local space for Fe
		linalg::types::Vector<Real, linalg::types::backend::CPU> Fe( (EvalEle::NodesPerElement * topoDOF.dofsPerNode()) );

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){
				
			// get element rngTags
			const Int* rngTags = mesh.getBoundaryTag(e);

			// face loop
			for (Index f = 0; f < mesh.data.facesPerElement; ++f){
				
				// get face rng tag
				Int rngTag = rngTags[f];
				const Index nodesPerFace = Basis::nodesPerFace(rngTag);

				Real faceNodeCoords[EvalEle::SpatialDim * nodesPerFace]
				Index faceNodeIDs[nodesPerFace];
				Basis::getFaceNodes(rngTag, faceNodeIDs);

				// extract face coordinates
				for (Index i = 0; i < nodesPerFace; ++i){
					
					const Real* faceNodeCoordsPtr = mesh.getNodeCoord(faceNodeIDs[i]);
					for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
						faceNodeCoords[EvalEle::SpatialDim*i + sD] = faceNodeCoordsPtr[sD];
					}
				
				}

				for (Index i = 0; i < nodesPerFace; ++i){
					for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){

						Index TdofIDi = topoDOF.getNodeDOF(faceNodeIDs[i], j);
						Index AdofIDi;
						if (topoDOF.isConstrained(TdofID) {
							continue;
						} else {
							AdofIDi = topoDOF.toAlgebraic(TdofIDi);
						}

						if (bcRegistry.isNatural(rngTag, j)){
							
							auto BoundaryConditions bcRegistry.getBCs(rngTag);

							for (auto& BoundaryCondition : BoundaryConditions){
								
								// get element data
								EvalEle evalE;
								evalE.bindElement(faceNodeCoords, time);
								
								// qp data
								EvalQP qp(EvalE);
								Real xi[Quadrature::NumPointsTotal*EvalEle::ParametricDim];
								Real w[Quadrature::NumPointsTotal];
								Quadrature::getPoints(xi);
								Quadrature::getWeights(w);

								// quadrature loop
								for (Index q = 0; q < Quadrature::NumPointsTotal; ++q){
									qp.evaluate(&xi[EvalEle::ParametricDim*q], w[q]);
									form.computeElementLevelVector(qp, Fe.data());
								}

								// scatter Fe into F
								for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
									for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
										
										Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
										if (topoDOF.isConstrained(TdofIDi)) continue;
										Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
										
										F.data()[AdofIDi] += Fe.data()[i*(topoDOF.dofsPerNode()) + j];

									}
								}

							} else {
								continue;
							}
						
						}

					}
				}

			}

		}

	}

} // namespace pdesolver::fem::boundary
