#ifndef POISSON_EVALQUADRATUREPOINTBOUNDARY_HPP
#define POISSON_EVALQUADRATUREPOINTBOUNDARY_HPP

#include "fem/eval/EvalQuadraturePointBoundary.hpp"

namespace pdesolver::fem::eval {

	template<typename Element, typename Basis, typename Geometry>
	class PoissonEvalQuadraturePointBoundary {
	public:

		Element element;
		Int faceID;

		PoissonEvalQuadraturePointBoundary(const Element& elem, const Int fID) : element(elem), faceID(fID) {}

		// dimensions
		static constexpr Index NodesPerElement = Element::NodesPerElement;
		static constexpr Int SpatialDim = Element::SpatialDim;
		static constexpr Int ParametricDim = Element::ParametricDim;

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

		// physical coordinate
		Real x[SpatialDim];

		// quadrature
		Real xi[ParametricDim-1];
		Real w;

		// ref basis values
		Real N[Element::NodesPerElement];

		// ref gradients
		Real dNdxi[ParametricDim*Element::NodesPerElement];

		// normal vector
		Real normal[SpatialDim];

		// geometry
		Real J[SpatialDim*ParametricDim];
		Index tangentID[ParametricDim-1];
		Real nCoeff;

		PDE_HOST PDE_DEVICE void evaluate(const Real* xi_q, const Real weight){
			
			// set quad info (need to add missing coordinate here, use faceID info)
			for (Index pD = 0; pD < ParametricDim; ++pD){
				if (){
					xi[pD] = 
				} else {
					xi[pD] = xi_q[pD];
				}
			}
			w = weight;
			
			// evaluate basis
			Basis::eval(xi, N);
			Basis::evalGradient(xi, dNdxi);
			nCoeff = Basis::getFaceTopology(faceID, tangentID);

			// geometry
			Geometry::mapToPhysical(coords, N, x);
			Geometry::computeJacobian(coords, dNdxi, J);
			Geometry::computeNormal(J, tangentID, nCoeff, normal);

		}

	}; // class PoissonEvalQuadraturePointBoundary

} // namespace pdesolver::fem::eval

#endif
