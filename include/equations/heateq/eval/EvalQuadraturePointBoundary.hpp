#ifndef HEATEQUATION_EVALQUADRATUREPOINTBOUNDARY_HPP
#define HEATEQUATION_EVALQUADRATUREPOINTBOUNDARY_HPP

#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/eval/EvalQuadraturePointBoundary.hpp"

namespace pdesolver::equations::heateq {

	template<typename Element, typename Basis, typename Geometry>
	class EvalQuadraturePointBoundary {
	public:

		static constexpr Index SpatialDim = Element::SpatialDim;
		static constexpr Index ParametricDim = Element::ParametricDim;

		Element element;
		Int faceID;

		Index faceNodeLocalIDs[fem::dispatch::kMaxNodesPerElementBoundary<ParametricDim>];

		EvalQuadraturePointBoundary(const Element& elem, const Int fID) : element(elem), faceID(fID) {

			element.basis().getFaceNodes(faceID, faceNodeLocalIDs);

		}

		Index nodesPerElement() const { return element.nodesPerElement(); }

		Index nodesPerFace() const { return element.basis().nodesPerFace(faceID); }

		// parent element attributes
		const Real time = element.t;
		const Real* coords = element.nodeCoords;

		// physical coordinates
		Real x[SpatialDim];

		// reference coordinate
		Real xi[ParametricDim];

		// quadrature
		Real xi_face[ParametricDim-1];
		Real w;

		Real N[fem::dispatch::kMaxNodesPerElement<ParametricDim>];

		Real Nface[fem::dispatch::kMaxNodesPerElementBoundary<ParametricDim>];

		// ref gradients
		Real dNdxi[ParametricDim*fem::dispatch::kMaxNodesPerElement<ParametricDim>];

		// normal vectors
		Real normal[SpatialDim];
		Real normalRef[ParametricDim];

		// geometry
		Real J[SpatialDim*ParametricDim];

		PDE_HOST PDE_DEVICE void evaluate(const Real* xi_face_q, const Real weight){

			// set quad info
			for (Index pD = 0; pD < (ParametricDim - 1); ++pD){
				xi_face[pD] = xi_face_q[pD];
			}
			w = weight;

			// get volume parametric coordinates
			element.basis().mapFaceToElement(faceID, xi_face, xi);
			element.basis().eval(xi, N);
			element.basis().evalGradient(xi, dNdxi);
			element.basis().getFaceTopology(faceID, normalRef);

			// gather N onto face-local indices
			for (Index a = 0; a < nodesPerFace(); ++a){
				Nface[a] = N[faceNodeLocalIDs[a]];
			}

			// geometry
			Geometry::mapToPhysical(coords, N, x, nodesPerElement());
			Geometry::computeJacobian(coords, dNdxi, J, nodesPerElement());
			Geometry::computeBoundaryNormal(J, normalRef, normal);

		}

	}; // class EvalQuadraturePointBoundary

} // namespace pdesolver::equations::heateq

#endif
