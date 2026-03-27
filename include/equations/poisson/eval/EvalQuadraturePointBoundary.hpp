#ifndef POISSON_EVALQUADRATUREPOINTBOUNDARY_HPP
#define POISSON_EVALQUADRATUREPOINTBOUNDARY_HPP

#include "fem/eval/EvalQuadraturePointBoundary.hpp"

namespace pdesolver::fem::eval {

	template<typename Element, typename Basis, typename Geometry>
	class PoissonEvalQuadraturePointBoundary {
	public:

		Element element;
		Int faceID;

		PoissonEvalQuadraturePointBoundary(const Element& elem, const Int rngID) : element(elem), faceID(rngID) {}

		// dimensions
		static constexpr Index NodesPerFace = Element::NodesPerFace;
		static constexpr Int SpatialDim = Element::SpatialDim;
		static constexpr Int ParametricDim = Element::ParametricDim;
	
		// parent element attributes
		const Real time = element.t;
		const Real* coords = element.nodeCoords;

		// physical coordinate
		Real x[SpatialDim];

		// quadrature
		Real xi[ParametricDim];
		Real w;

		// ref basis values
		Real N[NodesPerFace];

		// ref gradients
		Real dNdxi[ParametricDim*NodesPerFace];

		// normal vector
		Real normal[SpatialDim];

		// geometry
		Real J[SpatialDim*ParametricDim];
		Index tangentID[ParametricDim - 1];
		Real nCoeff;

		PDE_HOST PDE_DEVICE void evaluate(const Real* xi_q, const Real weight){
			
			// set quad info
			for (Index pD = 0; pD < ParametricDim; ++pD){
				xi[pD] = xi_q[pD];
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
