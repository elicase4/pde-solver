namespace residuum::fem::assembly {

	template<Index numDOFs, Index SpatialDim, GatherMode Mode, typename VectorT>
	PDE_HOST PDE_DEVICE void gatherElementVector(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry* bcRegistry, const Real time, Real* Ue, const VectorT* U){

		if constexpr (Mode == GatherMode::Free) {

			for (Index i = 0; i < nodesPerElement; ++i){
				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

					Ue[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j] = U->data()[AdofIDi];

				}
			}

		} else {

			for (Index i = 0; i < nodesPerElement; ++i){

				Real bcVal[topology::TopologicalDOF<numDOFs>::dofsPerNode];
				bool haveBcVal = false;

				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					
					// gather free nodes
					if (!topoDOF.isConstrained(TdofIDi)) {
						if constexpr (Mode == GatherMode::Full) {
							Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
							Ue[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j] = U->data()[AdofIDi];
						}
						continue;
					}

					// gather constrained nodes
					if (!haveBcVal) {
						Int rngTag = topoDOF.getConstraintTag(TdofIDi);
						const auto* entries = bcRegistry->getEntries(rngTag);
						if (entries) {
							for (const auto& entry : *entries) {
								entry->eval(time, &nodeCoords[SpatialDim*i], bcVal);
							}
						}
						haveBcVal = true;
					}

					Ue[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j] = bcVal[j];

				}

			}

		}

	}

	template<Index numDOFs, ScatterMode Mode, typename VectorT>
	PDE_HOST PDE_DEVICE void scatterElementVector(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Xe, VectorT& X, Real coefficient){

		for (Index i = 0; i < nodesPerElement; ++i){
			for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
				if (topoDOF.isConstrained(TdofIDi)) continue;
				Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

				X.data()[AdofIDi] += coefficient * Xe[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j];

			}
		}

	}

	template<Index numDOFs, ScatterMode Mode, typename MatrixT>
	PDE_HOST PDE_DEVICE void scatterElementMatrix(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Ke, MatrixT& K, Real coefficient){

		for (Index i = 0; i < nodesPerElement; ++i){
			for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
				if (topoDOF.isConstrained(TdofIDi)) continue;
				Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

				for (Index k = 0; k < nodesPerElement; ++k){
					for (Index l = 0; l < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++l){

						Index TdofIDk = topoDOF.getNodeDOF(nodeIDs[k], l);
						if (topoDOF.isConstrained(TdofIDk)) continue;
						Index AdofIDk = topoDOF.toAlgebraic(TdofIDk);
						Index p = K.getDataIndex(AdofIDi, AdofIDk);

						K.data()[p] += coefficient * Ke[(i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j)*(nodesPerElement * topology::TopologicalDOF<numDOFs>::dofsPerNode) + (k*topology::TopologicalDOF<numDOFs>::dofsPerNode + l)];

					}
				}

			}
		}

	}

} // namespace residuum::fem::assembly
