namespace pdesolver::fem::boundary {

template<>
class BoundaryApplicator<linalg::types::backend::CPU> {
public:
	
	template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature>
	void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const EssentialBoundaryRegistry& bcRegistry, const Real time, const Model& model, const FormRegistry& forms, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate Fe on the stack
		Real Fe[(EvalEle::NodesPerElement*topology::TopologicalDOF<numDOFs>::dofsPerNode)];

		// allocate Ge on the stack
		Real Ge[(EvalEle::NodesPerElement*topology::TopologicalDOF<numDOFs>::dofsPerNode)];

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){
		
			// zero-out Fe
			std::memset(Fe, 0.0, sizeof(Fe));

			// zero-out Ge
			std::memset(Ge, 0.0, sizeof(Ge));
			
			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * EvalEle::NodesPerElement];

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

				Real bcVal[topology::TopologicalDOF<numDOFs>::dofsPerNode];

				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){
					
					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (!topoDOF.isConstrained(TdofIDi)) {
						continue;
					}
					
					Int rngTag = topoDOF.getConstraintTag(TdofIDi);

					const auto* entries = bcRegistry.getEntries(rngTag);

					if (entries) {
						for (const auto& entry : *entries) {
							entry->eval(time, &nodeCoords[EvalEle::SpatialDim*i], bcVal);
						}
					}

					Ge[i*topology::TopologicalDOF<numDOFs>::dofsPerNode+j] = bcVal[j];

				}

			}


			// quadrature loop
			for (Index q = 0; q < Quadrature::NumPointsTotal; ++q){
				qp.evaluate(&xi[EvalEle::ParametricDim*q], w[q]);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelVector(qp, Ge, Fe);
			}

			// scatter Fe into F
			for (Index i = 0; i < EvalEle::NodesPerElement; ++i){
				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){
					
					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
					
					F.data()[AdofIDi] -= Fe[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j];

				}
			}

		}

	}

	template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Quadrature>
	void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const NaturalBoundaryRegistry<EvalQP>& bcRegistry, const Real time, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate local space for Fe
		Real Fe[(EvalEle::NodesPerElement*topology::TopologicalDOF<numDOFs>::dofsPerNode)];

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){
				
			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim*EvalEle::NodesPerElement];

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
				std::memset(Fe, 0.0, sizeof(Fe));

				// get face rng tag
				Int rngTag = rngTags[f];
				if (rngTag < 0) continue;
				if (!bcRegistry.hasAny(rngTag)) continue;
				const Index nodesPerFace = EvalQP::NodesPerFace(rngTag);

				Real faceNodeCoords[EvalEle::SpatialDim*EvalQP::NodesPerElement] = {0};
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
				Real bcVal[numDOFs];

				// quadrature loop
				for (Index q = 0; q < Quadrature::NumPointsTotal; ++q){
					
					qp.evaluate(&xi[(EvalEle::ParametricDim-1)*q], w[q]);
					bcRegistry.apply(rngTag, qp, time, faceNodeCoords, bcVal, Fe);
				
				}

				// scatter Fe into F
				for (Index i = 0; i < nodesPerFace; ++i){
					
					for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){
						
						Index TdofIDi = topoDOF.getNodeDOF(elemNodeGlobalIDs[faceNodeLocalIDs[i]], j);
						if (topoDOF.isConstrained(TdofIDi)) continue;
						Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
						
						F.data()[AdofIDi] += Fe[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j];

					}
				
				}
			
			}
		
		}

		
	}

}; // class BoundaryApplicator

} // namespace pdesolver::fem::boundary
