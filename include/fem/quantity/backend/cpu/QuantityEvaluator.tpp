#include <cmath>
#include <cstring>

#include "fem/assembly/Assembler.hpp"

namespace pdesolver::fem::quantity {

template<>
class QuantityEvaluator<linalg::types::backend::CPU> {
public:

	template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename QuantityFormsT, typename Quadrature>
	static void evaluateDomain(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const Model& model, const QuantityFormsT& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, const std::array<const linalg::types::Vector<Real, linalg::types::backend::CPU>*, EvalQP::NumAuxStates>& auxStates, Real* out){

		// zero-out output buffer & measure summation
		for (Index c = 0; c < QuantityFormsT::TotalComponents; ++c) out[c] = 0.0;
		Real measureSum = 0.0;

		// allocate Ue on the stack
		Real Ue[fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim>*numDOFs];

		// allocate the per-element auxiliary-state buffers on the stack
		Real Ue_aux[EvalQP::NumAuxStates][fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim>*numDOFs];
		const Real* Ue_auxPtrs[EvalQP::NumAuxStates];

		EvalEle localEle = evalEle;

		// quadrature points/weights
		Real xi[fem::dispatch::kMaxQuadraturePointsTotal<EvalEle::ParametricDim>*EvalEle::ParametricDim];
		Real w[fem::dispatch::kMaxQuadraturePointsTotal<EvalEle::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e) {

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim * fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			// gather the complete nodal solution
			fem::assembly::gatherElementVector<numDOFs, EvalEle::SpatialDim, fem::assembly::GatherMode::Full>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ue, &U);

			// gather each auxiliary vector into its own row -- always Free: a rate/derivative field
			// has nothing meaningful to look up via the essential BC's VALUE eval() function at a
			// constrained node (matches Assembler::assembleMatrix/assembleVector's aux gather)
			for (Index s = 0; s < EvalQP::NumAuxStates; ++s) {
				std::memset(Ue_aux[s], 0.0, sizeof(Ue_aux[s]));
				if (auxStates[s] != nullptr) {
					fem::assembly::gatherElementVector<numDOFs, EvalEle::SpatialDim, fem::assembly::GatherMode::Free>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ue_aux[s], auxStates[s]);
				}
				Ue_auxPtrs[s] = Ue_aux[s];
			}

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// qp data
			EvalQP qp(localEle);

			// quadrature loop
			for (Index q = 0; q < quadrature.numPointsTotal(); ++q){

				qp.evaluate(&xi[EvalEle::ParametricDim*q], w[q]);
				qp.interpolateFields(Ue, Ue_auxPtrs);
				model.eval(qp);
				model.evalGradient(qp);

				forms.computeElementLevelValue(qp, Ue, out);

				measureSum += qp.measure * qp.w;

			}

		}

		// finalize Average-mode slices
		forms.finalize(out, measureSum);

	}

	template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename QuantityFormsT, typename Quadrature>
	static void evaluateBoundaryRegistry(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const Model& model, const QuantityFormsT& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, linalg::types::backend::CPU>& U, BoundaryQuantityRegistry<QuantityFormsT>& registry){

		// zero every registered tag's slot & measureSum
		for (auto& [tag, entry] : registry.entries()) {
			entry.result.fill(0.0);
			entry.measureSum = 0.0;
		}

		// allocate Ue on the stack
		Real Ue[fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim>*numDOFs];

		EvalEle localEle = evalEle;

		// boundary quadrature
		Real xi[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEle::ParametricDim>*(EvalEle::ParametricDim-1)];
		Real w[fem::dispatch::kMaxQuadraturePointsTotalBoundary<EvalEle::ParametricDim>];
		quadrature.getPoints(xi);
		quadrature.getWeights(w);

		// element loop
		for (Index e = 0; e < mesh.data.numElements; ++e) {

			// get element rngTags, skip elements with no face on ANY registered tag at all
			const Int* rngTags = mesh.getBoundaryTag(e);

			bool touchesAny = false;
			for (Index f = 0; f < mesh.data.facesPerElement; ++f) {
				if (registry.hasTag(rngTags[f])) { touchesAny = true; break; }
			}
			if (!touchesAny) continue;

			// extract node coordinates
			const Index* nodeIDs = mesh.getElementNodes(e);
			Real nodeCoords[EvalEle::SpatialDim*fem::dispatch::kMaxNodesPerElement<EvalEle::ParametricDim>];

			for (Index i = 0; i < localEle.nodesPerElement(); ++i){

				const Real* nodeCoordsPtr = mesh.getNodeCoord(nodeIDs[i]);

				for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD){
					nodeCoords[EvalEle::SpatialDim*i + sD] = nodeCoordsPtr[sD];
				}

			}

			fem::assembly::gatherElementVector<numDOFs, EvalEle::SpatialDim, fem::assembly::GatherMode::Full>(nodeIDs, localEle.nodesPerElement(), nodeCoords, topoDOF, &bcRegistry, time, Ue, &U);

			// bind element data
			localEle.bindElement(nodeCoords, time);

			// face loop -- each face accumulates into ITS OWN tag's slot, if registered
			for (Index f = 0; f < mesh.data.facesPerElement; ++f) {

				if (!registry.hasTag(rngTags[f])) continue;

				auto& entry = registry.entries().at(rngTags[f]);

				// qp data
				EvalQP qp(localEle, f);

				// quadrature loop
				for (Index q = 0; q < quadrature.numPointsTotal(); ++q){

					qp.evaluate(&xi[(EvalEle::ParametricDim-1)*q], w[q]);
					qp.interpolateFields(Ue);
					model.eval(qp);
					model.evalGradient(qp);

					forms.computeElementLevelValue(qp, Ue, entry.result.data());

					Real normalMag = 0;
					for (Index sD = 0; sD < EvalEle::SpatialDim; ++sD) {
						normalMag += qp.normal[sD] * qp.normal[sD];
					}
					entry.measureSum += std::sqrt(normalMag) * qp.w;

				}

			}

		}

		// finalize Average-mode slices per tag, using each tag's OWN measureSum
		for (auto& [tag, entry] : registry.entries()) {
			forms.finalize(entry.result.data(), entry.measureSum);
		}

	}

}; // class QuantityEvaluator<linalg::types::backend::CPU>

} // namespace pdesolver::fem::quantity
