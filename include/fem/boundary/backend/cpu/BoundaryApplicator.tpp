namespace pdesolver::fem::boundary {

template<>
class BoundaryApplicator<linalg::types::backend::CPU> {
public:
	
	template<eval::EvalElement EvalEle, typename EvalQP, typename Form, typename Model, typename Quadrature, typename Function>
	static void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& bcRegistry, const Real time, const Model& model, const Form& form, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate local space for Fe
		linalg::types::Vector<Real, linalg::types::backend::CPU> Fe( (EvalEle::NodesPerElement * topoDOF.dofsPerNode()) );

		// allocate local space for Ge
		linalg::types::Vector<Real, linalg::types::backend::CPU> Ge( (EvalEle::NodesPerElement * topoDOF.dofsPerNode()) );

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){
		
			// zero-out Fe
			Fe.zero();
			
			// zero-out Ge
			Ge.zero();
			
			// extract node coordinates
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
			EvalQP qp(evalE);
			Real xi[Quadrature::NumPointsTotal*EvalEle::ParametricDim];
			Real w[Quadrature::NumPointsTotal];
			Quadrature::getPoints(xi);
			Quadrature::getWeights(w);
			
			// fill Ge
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				
				Real bcVal[BoundaryCondition<Function>::NumComponents];
				Int rngTag;
				
				for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) {
						rngTag = topoDOF.getConstraintTag(TdofIDi);
					} else {
						continue;
					}
					
					// extract and evaluate boundary value function
					for (const auto& entry: bcRegistry.entries()){
						if (entry->tag() != rngTag) continue;
						if (entry->componentType(j) != BCCategory::Essential) continue;
						entry->eval(qp.time, qp.x, bcVal);
					}
					
					if (bcRegistry.isEssential(rngTag, j)){
						
						Ge.data()[i*topoDOF.dofsPerNode() + j] = bcVal[j];
					
					}

				}

			}

			// quadrature loop
			for (Index q = 0; q < Quadrature::NumPointsTotal; ++q){
				qp.evaluate(&xi[EvalEle::ParametricDim*q], w[q]);
				model.eval(qp);
				model.evalGradient(qp);
				form.computeElementLevelVector(qp, Ge.data(), Fe.data());
			}
			
			// scatter Fe into F
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
					
					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
					
					F.data()[AdofIDi] -= Fe.data()[i*topoDOF.dofsPerNode() + j];

				}
			}

		}

	}

	template<eval::EvalElement EvalEle, typename EvalQP, typename Form, typename Quadrature>
	static void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& bcRegistry, const Real time, const Form& form, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate local space for Fe
		linalg::types::Vector<Real, linalg::types::backend::CPU> Fe( (EvalEle::NodesPerElement * topoDOF.dofsPerNode()) );

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){
				
			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * EvalEle::NodesPerElement];

			// extract node coordinates
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}
			
			// get element rngTags
			const Int* rngTags = mesh.getBoundaryTag(e);

			// face loop
			for (Index f = 0; f < mesh.data.facesPerElement; ++f){
				
				// zero-out Fe
				Fe.zero();

				// get face rng tag
				Int rngTag = rngTags[f];
				if (rngTag < 0) continue;
				if (!bcRegistry.hasAny(rngTag)) continue;
				const Index nodesPerFace = EvalQP::NodesPerFace(rngTag);

				Real faceNodeCoords[EvalEle::SpatialDim * EvalQP::NodesPerElement];
				Index faceNodeLocalIDs[EvalQP::NodesPerElement];
				EvalQP::getFaceNodes(rngTag, faceNodeLocalIDs);
				const Index* elemNodeGlobalIDs = mesh.getElementNodes(e);

				// extract face coordinates
				for (Index i = 0; i < nodesPerFace; ++i){
					
					Index faceNodeGlobalID = elemNodeGlobalIDs[faceNodeLocalIDs[i]];
					const Real* faceNodeCoordsPtr = mesh.getNodeCoord(faceNodeGlobalID);
					
					for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
						faceNodeCoords[EvalEle::SpatialDim*i + sD] = faceNodeCoordsPtr[sD];
					}
				
				}

				// get element data
				EvalEle evalE;
				evalE.bindElement(nodeCoords, time);
				
				// qp data
				EvalQP qp(evalE, f, faceNodeCoords);
				Real xi[Quadrature::NumPointsTotal*(EvalEle::ParametricDim-1)];
				Real w[Quadrature::NumPointsTotal];
				Quadrature::getPoints(xi);
				Quadrature::getWeights(w);

				// quadrature loop
				for (Index q = 0; q < Quadrature::NumPointsTotal; ++q){
					qp.evaluate(&xi[(EvalEle::ParametricDim-1)*q], w[q]);
					form.computeElementLevelVector(qp, nullptr, Fe.data());
				}

				// scatter Fe into F
				for (Index i = 0; i < nodesPerFace; ++i){
					
					for (Index j = 0; j < topoDOF.dofsPerNode(); ++j){
						
						if (!bcRegistry.isNatural(rngTag, j)) continue;
						Index TdofIDi = topoDOF.getNodeDOF(elemNodeGlobalIDs[faceNodeLocalIDs[i]], j);
						if (topoDOF.isConstrained(TdofIDi)) continue;
						Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
						
						F.data()[AdofIDi] += Fe.data()[i*topoDOF.dofsPerNode() + j];

					}
				}
			
			}
		
		}

		
	}

}; // class BoundaryApplicator

} // namespace pdesolver::fem::boundary
