#ifndef HEATEQUATION_EVALQUADRATUREPOINTVOLUME_HPP
#define HEATEQUATION_EVALQUADRATUREPOINTVOLUME_HPP

#include "fem/DiscretizationLimits.hpp"
#include "fem/eval/EvalQuadraturePointVolume.hpp"

namespace pdesolver::equations::heateq {

	template<typename Element, typename Basis, typename Geometry>
	class EvalQuadraturePointVolume {
	public:

		Element element;

		EvalQuadraturePointVolume(const Element& elem) : element(elem) {}

		static constexpr Index SpatialDim = Element::SpatialDim;
		static constexpr Index ParametricDim = Element::ParametricDim;

		Index nodesPerElement() const { return element.nodesPerElement(); }

		// parent element attributes
		const Real time = element.t;
		const Real* coords = element.nodeCoords;

		// physical coordinate
		Real x[SpatialDim];

		// quadrature
		Real xi[ParametricDim];
		Real w;

		Real N[fem::kMaxNodesPerElement<ParametricDim>];

		// ref gradients
		Real dNdxi[ParametricDim*fem::kMaxNodesPerElement<ParametricDim>];

		// physical gradients
		Real dNdx[SpatialDim*fem::kMaxNodesPerElement<ParametricDim>];

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

			element.basis().eval(xi, N);
			element.basis().evalGradient(xi, dNdxi);

			Geometry::mapToPhysical(coords, N, x, nodesPerElement());
			Geometry::computeJacobian(coords, dNdxi, J, nodesPerElement());
			Geometry::computeMetric(J, g);
			measure = Geometry::computeMeasure(g);

			// transforms
			Geometry::transformGradient(J, g, dNdxi, dNdx, nodesPerElement());

		}

	}; // class EvalQuadraturePointVolume

} // namespace pdesolver::equations::heateq

#endif
