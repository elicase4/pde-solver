namespace pdesolver::fem::boundary {

template<>
class BoundaryApplicator<linalg::types::backend::CPU>
public:
	
	template<eval::EvalElement EvalEle, typename EvalQP, typename Form, typename Basis, typename Model>
	void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& bcRegistry, const Real time, const Model& model, const Form&form, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate local space for Fe
		linalg::types::Vector<Real, linalg::types::backend::CPU> Fe( (EvalEle::NodesPerElement * topoDOF.dofsPerNode()) );

		// allocate local space for Ge
		linalg::types::Vector<Real, linalg::types::backend::CPU> Ge( (EvalEle::NodesPerElement * topoDOF.dofsPerNode()) );

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){
		
			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * EvalEle::NodesPerElement];

			// extract node coordinates and fill Ge
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// zero-out Fe
			Fe.zero();
			
			// zero-out Ge
			Ge.zero();
			
			// gather G into Ge
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				// fix this loop to compensate for bcs with multiple dof components
				for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					Int rngTag;
					if (topoDOF.isConstrained(TdofIDi){
						rngTag = topoDOF.getConstraintTag(TdofIDi);
					} else {
						continue;
					}
					
					if (bcRegistry.isEssential(rngTag, j)){
						
						auto BoundaryConditions bcRegistry.getBCs(rngTag);

						for (auto& BoundaryCondition : BoundaryConditions){
				
							Real outValue[BoundaryCondition.f.NumComponents()];
							BoundaryCondition.f.eval(time, &nodeCoords[i*EvalEle::NodesPerElement], outValue);
							
							for (k = 0; k < BoundaryCondition.f.NumComponents(); ++k){
								Ge.data()[i*(EvalElel::NodesPErElement*topoDOF.dofsPerNode()) + j + k] = outValue[k];
							}

						}
					
					}

				}
			}

		}

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

				// zero-out Fe
				Fe.zero();

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
								for (Index i = 0; i < NodesPerFace; ++i){
									for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
										
										Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
										if (topoDOF.isConstrained(TdofIDi)) continue;
										Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
										
										F.data()[AdofIDi] += Fe.data()[i*(nodesPerFace*topoDOF.dofsPerNode()) + j];

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
