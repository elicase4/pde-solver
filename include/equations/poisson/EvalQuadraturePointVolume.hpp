#ifndef POISSON_EVALQUADRATUREPOINTVOLUME_HPP
#define POISSON_EVALQUADRATUREPOINTVOLUME_HPP

#include "fem/eval/EvalQuadraturePointVolume.hpp"

namespace pdesolver::fem::eval {

	template<typename Element, typename Basis, typename Geometry>
	class PoissonEvalQuadraturePointVolume {
	public:

		Element element;

		PoissonEvalQuadraturePointVolume(const Element& elem) : element(elem) {}

		// dimensions
		static constexpr Index NodesPerElement = Element::NodesPerElement;
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
		Real N[NodesPerElement];

		// ref gradients
		Real dNdxi[ParametricDim*NodesPerElement];

		// physical gradients
		Real dNdx[SpatialDim*NodesPerElement];

		// geometry
		Real J[SpatialDim*ParametricDim];
		Real g[ParametricDim*ParametricDim];

		// measure
		Real measure;

		// conductivity coefficient
		Real K[SpatialDim*SpatialDim];

		PDE_HOST PDE_DEVICE void evaluate(const Real* xi_q, const Real weight){
			
			// set quad info
			for (Index pD = 0; pD < ParametricDim; ++pD){
				xi[pD] = xi_q[pD];
			}
			w = weight;
			
			// evaluate basis
			Basis::eval(xi, N);
			Basis::evalGradient(xi, dNdxi);

			// geometry
			Geometry::mapToPhysical(coords, N, x);
			Geometry::computeJacobian(coords, dNdxi, J);
			Geometry::computeMetric(J, g);
			measure = Geometry::computeMeasure(g);

			// transforms
			Geometry::transformGradient(J, g, dNdxi, dNdx);

		}

	}; // class PoissonEvalQuadraturePointVolume

} // namespace pdesolver::fem::eval

#endif
