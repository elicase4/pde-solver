namespace pdesolver::fem::boundary {

template<>
class BoundaryApplicator<linalg::types::backend::CPU> {
public:

	template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature>
	void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const EssentialBoundaryRegistry& bcRegistry, const Real time, const Model& model, const FormRegistry& forms, const EvalEle& evalEle, const Quadrature& quadrature, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate Fe on the stack
		Real Fe[(fem::kMaxNodesPerElement<EvalEle::ParametricDim>*topology::TopologicalDOF<numDOFs>::dofsPerNode)];

		// allocate Ge on the stack
		Real Ge[(fem::kMaxNodesPerElement<EvalEle::ParametricDim>*topology::TopologicalDOF<numDOFs>::dofsPerNode)];

		EvalEle localEle = evalEle;

		// quadrature points/weights
		Real xi[fem::kMaxQuadraturePointsTotal<EvalEle::ParametricDim>*EvalEle::ParametricDim];
		Real w[fem::kMaxQuadraturePointsTotal<EvalEle::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){

			// zero-out Fe
			std::memset(Fe, 0.0, sizeof(Fe));

			// zero-out Ge
			std::memset(Ge, 0.0, sizeof(Ge));

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * fem::kMaxNodesPerElement<EvalEle::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQP qp(localEle);

			// fill Ge
			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

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
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEle::ParametricDim*q], w[q]);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelVector(qp, Ge, Fe);
			}

			// scatter Fe into F
			for (Index i = 0; i < localEle.nodesPerElement(); ++i){
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
	static void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const NaturalBoundaryRegistry<EvalQP>& bcRegistry, const Real time, const EvalEle& evalEle, const Quadrature& quadrature, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate local space for Fe
		Real Fe[(fem::kMaxNodesPerElement<EvalEle::ParametricDim>*topology::TopologicalDOF<numDOFs>::dofsPerNode)];

		EvalEle localEle = evalEle;

		// boundary quadrature
		Real xi[fem::kMaxQuadraturePointsTotalBoundary<EvalEle::ParametricDim>*(EvalEle::ParametricDim-1)];
		Real w[fem::kMaxQuadraturePointsTotalBoundary<EvalEle::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim*fem::kMaxNodesPerElement<EvalEle::ParametricDim>];

			// extract node coordinates
			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// bind element data
			localEle.bindElement(nodeCoords, time);

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

				// get face nodes information
				const Index nodesPerFace = localEle.basis().nodesPerFace(f);
				
				Real faceNodeCoords[EvalEle::SpatialDim*fem::kMaxNodesPerElement<EvalEle::ParametricDim>] = {0};
				Index faceNodeLocalIDs[fem::kMaxNodesPerElement<EvalEle::ParametricDim>];
				
				localEle.basis().getFaceNodes(f, faceNodeLocalIDs);
				const Index* elemNodeGlobalIDs = mesh.getElementNodes(e);

				// extract face coordinates
				for (Index i = 0; i < nodesPerFace; ++i){

					Index faceNodeGlobalID = elemNodeGlobalIDs[faceNodeLocalIDs[i]];
					const Real* faceNodeCoordsPtr = mesh.getNodeCoord(faceNodeGlobalID);

					for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
						faceNodeCoords[EvalEle::SpatialDim*i + sD] = faceNodeCoordsPtr[sD];
					}

				}

				// qp data
				EvalQP qp(localEle, f);
				Real bcVal[numDOFs*EvalEle::SpatialDim];

				// quadrature loop
				for (Index q = 0; q < quadrature.numPointsTotal(); ++q){

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
