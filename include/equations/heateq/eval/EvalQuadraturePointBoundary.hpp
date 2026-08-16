#ifndef HEATEQUATION_EVALQUADRATUREPOINTBOUNDARY_HPP
#define HEATEQUATION_EVALQUADRATUREPOINTBOUNDARY_HPP

#include "fem/DiscretizationLimits.hpp"
#include "fem/eval/EvalQuadraturePointBoundary.hpp"

namespace pdesolver::equations::heateq {

	template<typename Element, typename Basis, typename Geometry>
	class EvalQuadraturePointBoundary {
	public:

		static constexpr Index SpatialDim = Element::SpatialDim;
		static constexpr Index ParametricDim = Element::ParametricDim;

		Element element;
		Int faceID;

		Real faceCoords[SpatialDim*fem::kMaxNodesPerElement<ParametricDim>];

		Index faceNodeLocalIDs[fem::kMaxNodesPerElement<ParametricDim>];

		EvalQuadraturePointBoundary(const Element& elem, const Int fID, const Real* faceNodeCoords) : element(elem), faceID(fID) {

			element.basis().getFaceNodes(faceID, faceNodeLocalIDs);

			// set face coordinates
			for (Index a = 0; a < element.basis().nodesPerFace(faceID); ++a){
				for (Index sD = 0; sD < SpatialDim; ++sD){
					faceCoords[a*SpatialDim + sD] = faceNodeCoords[a*SpatialDim + sD];
				}
			}

		}

		Index nodesPerElement() const { return element.nodesPerElement(); }

		Index nodesPerFace() const { return element.basis().nodesPerFace(faceID); }

		// parent element attributes
		const Real time = element.t;
		const Real* coords = element.nodeCoords;

		// physical coordinates
		Real x[SpatialDim];
		Real x_face[SpatialDim];

		// reference coordinate
		Real xi[ParametricDim];

		// quadrature
		Real xi_face[ParametricDim-1];
		Real w;

		// ref basis values -- capped, see faceCoords above.
		Real N[fem::kMaxNodesPerElement<ParametricDim>];

		// ref gradients
		Real dNdxi[ParametricDim*fem::kMaxNodesPerElement<ParametricDim>];

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

			// geometry
			Geometry::mapToPhysical(coords, N, x, nodesPerElement());
			Geometry::mapToPhysical(faceCoords, N, x_face, nodesPerElement());
			Geometry::computeJacobian(coords, dNdxi, J, nodesPerElement());
			Geometry::computeBoundaryNormal(J, normalRef, normal);

		}

	}; // class EvalQuadraturePointBoundary

} // namespace pdesolver::equations::heateq

#endif
