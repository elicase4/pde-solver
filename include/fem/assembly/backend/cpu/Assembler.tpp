namespace pdesolver::fem::assembly {

template<>
class Assembler<linalg::types::backend::CPU> {
public:

	template<Index numDOFs>
	static linalg::types::CSRMatrix<Real, linalg::types::backend::CPU> createMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF){

		// allocate matrix
		linalg::types::CSRMatrix<Real, linalg::types::backend::CPU> K(topoDOF.numFreeDOFs(), topoDOF.numFreeDOFs());

		// TODO: make flat
		// allocate adjacency list
		std::vector<std::vector<Index>> adjList(topoDOF.numFreeDOFs());

		// build adjacency list
		for (Index e = 0; e < mesh.data.numElements; ++e){

			const Index* nodeIDs = mesh.getElementNodes(e);

			for (Index i = 0; i < mesh.data.nodesPerElement; ++i){
				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					Index AdofIDi;
					if (topoDOF.isConstrained(TdofIDi)) {
						continue;
					} else {
						AdofIDi = topoDOF.toAlgebraic(TdofIDi);
					}

					for (Index k = 0; k < mesh.data.nodesPerElement; ++k){
						for (Index l = 0; l < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++l){

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

	template<Index numDOFs>
	static linalg::types::Vector<Real, linalg::types::backend::CPU> createVector(const mesh::Mesh&, const topology::TopologicalDOF<numDOFs>& topoDOF){

		linalg::types::Vector<Real, linalg::types::backend::CPU> F(topoDOF.numFreeDOFs());
		return F;

	}

	template<Index numDOFs, Index SpatialDim>
	static void gatherElementSolution(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, Real* Ue){

		for (Index i = 0; i < nodesPerElement; ++i){

			Real bcVal[topology::TopologicalDOF<numDOFs>::dofsPerNode];
			bool haveBcVal = false;

			for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

				Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);

				// free nodes
				if (!topoDOF.isConstrained(TdofIDi)) {
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);
					Ue[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j] = U.data()[AdofIDi];
					continue;
				}

				// constrained nodes
				if (!haveBcVal) {
					Int rngTag = topoDOF.getConstraintTag(TdofIDi);
					const auto* entries = bcRegistry.getEntries(rngTag);
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

	template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature>
	static void assembleMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const Model& model, const FormRegistry& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, linalg::types::CSRMatrix<Real, linalg::types::backend::CPU>& K){

		// allocate Ke on the stack
		Real Ke[(fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim> * numDOFs) * (fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim> * numDOFs)];

		// allocate Ue on the stack
		Real Ue[(fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim> * numDOFs)];

		// zero-out data in K
		K.zero();

		EvalEle localEle = evalEle;

		// setup quadrature points and weights
		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEle::ParametricDim>*EvalEle::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEle::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){

			// zero-out Ke
			std::memset(Ke, 0.0, sizeof(Ke));

			// zero-out Ue
			std::memset(Ue, 0.0, sizeof(Ue));

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// gather U into Ue
			for (Index i = 0; i < localEle.nodesPerElement(); ++i){
				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

					Ue[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j] = U.data()[AdofIDi];

				}
			}

			// gather any form-specific element data
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQP qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEle::ParametricDim*q], w[q]);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelMatrix(qp, Ue, Ke);
			}

			// scatter Ke into K
			for (Index i = 0; i < localEle.nodesPerElement(); ++i){
				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

					for (Index k = 0; k < localEle.nodesPerElement(); ++k){
						for (Index l = 0; l < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++l){

							Index TdofIDk = topoDOF.getNodeDOF(nodeIDs[k], l);
							if (topoDOF.isConstrained(TdofIDk)) continue;
							Index AdofIDk = topoDOF.toAlgebraic(TdofIDk);
							Index p = K.getDataIndex(AdofIDi, AdofIDk);

							K.data()[p] += Ke[(i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j)*(localEle.nodesPerElement() * topology::TopologicalDOF<numDOFs>::dofsPerNode) + (k*topology::TopologicalDOF<numDOFs>::dofsPerNode + l)];

						}
					}

				}
			}

		}

	}

	template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature>
	static void assembleVector(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const Model& model, const FormRegistry& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

		// allocate Fe on the stack
		Real Fe[fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim>*numDOFs];

		// allocate Ue on the stack
		Real Ue[fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim>*numDOFs];

		// zero-out data in F
		F.zero();

		// local mutable copy
		EvalEle localEle = evalEle;

		// quadrature points/weights
		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEle::ParametricDim>*EvalEle::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEle::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e){

			// zero-out Fe
			std::memset(Fe, 0.0, sizeof(Fe));

			// zero-out Ue
			std::memset(Ue, 0.0, sizeof(Ue));

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// gather U into Ue -- constrained DOFs left at 0 (no algebraic index for them in U).
			// Harmless today: nothing plugged into assembleVector (SourceForm/NodalSourceForm)
			// reads Ue at all. Would need Assembler::gatherElementSolution's complete gather
			// (constrained DOFs evaluated from an EssentialBoundaryRegistry) instead, the moment
			// a genuinely Ue-dependent LinearForm is wired up here (e.g. a Newton residual).
			for (Index i = 0; i < localEle.nodesPerElement(); ++i){
				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

					Ue[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j] = U.data()[AdofIDi];

				}
			}

			// gather any form-specific element data (no-op unless a registered form implements it)
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQP qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEle::ParametricDim*q], w[q]);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelVector(qp, Ue, Fe);
			}

			// scatter Fe into F
			for (Index i = 0; i < localEle.nodesPerElement(); ++i){
				for (Index j = 0; j < topology::TopologicalDOF<numDOFs>::dofsPerNode; ++j){

					Index TdofIDi = topoDOF.getNodeDOF(nodeIDs[i], j);
					if (topoDOF.isConstrained(TdofIDi)) continue;
					Index AdofIDi = topoDOF.toAlgebraic(TdofIDi);

					F.data()[AdofIDi] += Fe[i*topology::TopologicalDOF<numDOFs>::dofsPerNode + j];

				}
			}

		}

	}

}; // class Assembler <linalg::types::backend::CPU>

} // namespace pdesolver::fem::assembly
