namespace residuum::fem::assembly {

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

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, GatherMode Mode>
	requires evaluator::EvalModel<ModelT, EvalQPT>
	static void assembleMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, const std::array<const linalg::types::Vector<Real, linalg::types::backend::CPU>*, EvalQPT::NumAuxStates>& auxStates, linalg::types::CSRMatrix<Real, linalg::types::backend::CPU>& K, const fem::boundary::EssentialBoundaryRegistry* bcRegistry){

		// allocate Ke on the stack
		Real Ke[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs) * (fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];

		// allocate Ue on the stack
		Real Ue[(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];

		// allocate the per-element auxiliary-state buffers on the stack
		Real Ue_aux[EvalQPT::NumAuxStates][(fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim> * numDOFs)];
		const Real* Ue_auxPtrs[EvalQPT::NumAuxStates];

		// zero-out data in K
		K.zero();

		EvalEleT localEle = evalEle;

		// setup quadrature points and weights
		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*EvalEleT::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
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
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// gather U into Ue
			gatherElementVector<numDOFs, EvalEleT::SpatialDim, Mode>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue, &U);

			// aux states are rate fields, so a constrained node contributes 0 rather than a looked-up value
			for (Index s = 0; s < EvalQPT::NumAuxStates; ++s) {
				std::memset(Ue_aux[s], 0.0, sizeof(Ue_aux[s]));
				if (auxStates[s] != nullptr) {
					gatherElementVector<numDOFs, EvalEleT::SpatialDim, GatherMode::Free>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue_aux[s], auxStates[s]);
				}
				Ue_auxPtrs[s] = Ue_aux[s];
			}

			// gather any form-specific element data
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQPT qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				qp.interpolateFields(Ue, Ue_auxPtrs);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelMatrix(qp, Ke);
			}

			// scatter Ke into K
			scatterElementMatrix<numDOFs, ScatterMode::Free>(nodeIDs, localEle.nodesPerElement(), topoDOF, Ke, K);

		}

	}

	template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, GatherMode Mode>
	requires evaluator::EvalModel<ModelT, EvalQPT>
	static void assembleVector(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, const linalg::types::Vector<Real, linalg::types::backend::CPU>* fieldSource, const std::array<const linalg::types::Vector<Real, linalg::types::backend::CPU>*, EvalQPT::NumAuxStates>& auxStates, linalg::types::Vector<Real, linalg::types::backend::CPU>& F, const fem::boundary::EssentialBoundaryRegistry* bcRegistry){

		// allocate Fe on the stack
		Real Fe[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// allocate Ue on the stack
		Real Ue[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// allocate a buffer for fieldSource, used only when fieldSource != nullptr
		Real Ue_lin[fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];

		// allocate the per-element auxiliary-state buffers on the stack
		Real Ue_aux[EvalQPT::NumAuxStates][fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>*numDOFs];
		const Real* Ue_auxPtrs[EvalQPT::NumAuxStates];

		// zero-out data in F
		F.zero();

		// local mutable copy
		EvalEleT localEle = evalEle;

		// quadrature points/weights
		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>*EvalEleT::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEleT::ParametricDim>];
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
			Real nodeCoords[EvalEleT::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEleT::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEleT::SpatialDim; ++sD){
					nodeCoords[EvalEleT::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// gather U into Ue
			gatherElementVector<numDOFs, EvalEleT::SpatialDim, Mode>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue, &U);

			// gather fieldSource when it differs from the operand
			const Real* fieldPtr = Ue;
			if (fieldSource != nullptr) {
				std::memset(Ue_lin, 0.0, sizeof(Ue_lin));
				gatherElementVector<numDOFs, EvalEleT::SpatialDim, GatherMode::Full>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue_lin, fieldSource);
				fieldPtr = Ue_lin;
			}

			// aux states are rate fields, so a constrained node contributes 0 rather than a looked-up value
			for (Index s = 0; s < EvalQPT::NumAuxStates; ++s) {
				std::memset(Ue_aux[s], 0.0, sizeof(Ue_aux[s]));
				if (auxStates[s] != nullptr) {
					gatherElementVector<numDOFs, EvalEleT::SpatialDim, GatherMode::Free>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, bcRegistry, time, Ue_aux[s], auxStates[s]);
				}
				Ue_auxPtrs[s] = Ue_aux[s];
			}

			// gather any form-specific element data
			forms.gatherElementData(nodeIDs, localEle.nodesPerElement());

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQPT qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){
				qp.evaluate(&xi[EvalEleT::ParametricDim*q], w[q]);
				qp.interpolateFields(fieldPtr, Ue_auxPtrs);
				model.eval(qp);
				model.evalGradient(qp);
				forms.computeElementLevelVector(qp, Ue, Fe);
			}

			// scatter Fe into F
			scatterElementVector<numDOFs, ScatterMode::Free>(nodeIDs, localEle.nodesPerElement(), topoDOF, Fe, F);

		}

	}

}; // class Assembler <linalg::types::backend::CPU>

} // namespace residuum::fem::assembly
