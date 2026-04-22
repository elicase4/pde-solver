#ifndef POISSON_EVALQUADRATUREPOINTBOUNDARY_HPP
#define POISSON_EVALQUADRATUREPOINTBOUNDARY_HPP

#include "fem/eval/EvalQuadraturePointBoundary.hpp"

namespace pdesolver::fem::eval {

	template<typename Element, typename Basis, typename Geometry>
	class PoissonEvalQuadraturePointBoundary {
	public:

		// dimensions
		static constexpr Index NodesPerElement = Element::NodesPerElement;
		static constexpr Int SpatialDim = Element::SpatialDim;
		static constexpr Int ParametricDim = Element::ParametricDim;
		
		Element element;
		Int faceID;
		Real faceCoords[SpatialDim*Element::NodesPerElement];

		PoissonEvalQuadraturePointBoundary(const Element& elem, const Int fID, const Real* faceNodeCoords) : element(elem), faceID(fID) {
			
			// set face coordinates
			for (Index a = 0; a < NodesPerFace(faceID); ++a){
				for (Index sD = 0; sD < SpatialDim; ++sD){
					faceCoords[a*SpatialDim + sD] = faceNodeCoords[a*SpatialDim + sD];
				}
			}
			
		}

		// helper functions
		static Index NodesPerFace(const Int faceID){
			return Basis::nodesPerFace(faceID);
		}
	
		static void getFaceNodes(const Int faceID, Index* nodeIDs){
			Basis::getFaceNodes(faceID, nodeIDs);
		}
	
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

		// ref basis values
		Real N[Element::NodesPerElement];

		// ref gradients
		Real dNdxi[ParametricDim*Element::NodesPerElement];

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
			Basis::mapFaceToElement(faceID, xi_face, xi);
			
			// evaluate basis
			Basis::eval(xi, N);
			Basis::evalGradient(xi, dNdxi);
			Basis::getFaceTopology(faceID, normalRef);

			// geometry
			Geometry::mapToPhysical(coords, N, x);
			Geometry::mapToPhysical(faceCoords, N, x_face);
			Geometry::computeJacobian(coords, dNdxi, J);
			Geometry::computeBoundaryNormal(J, normalRef, normal);

		}

	}; // class PoissonEvalQuadraturePointBoundary

} // namespace pdesolver::fem::eval

#endif
