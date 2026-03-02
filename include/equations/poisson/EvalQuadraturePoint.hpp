#ifndef POISSON_EVALQUADRATUREPOINT_HPP
#define POISSON_EVALQUADRATUREPOINT_HPP

#include "fem/eval/EvalQuadraturePoint.hpp"

namespace pdesolver::fem::eval {

	template<typename Geometry, typename Basis>
	class PoissonEvalQuadraturePoint {
	public:

		// dimensions
		static constexpr Index NodesPerElement = Basis::NodesPerElement;
		static constexpr Int SpatialDim = Geometry::SpatialDim;
		static constexpr Int ParametricDim = Geometry::ParametricDim;
		
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

		// rhs function
		Real rhsF[SpatialDim];

		PDE_HOST PDE_DEVICE void evaluate(const Real* coords, const Real* xi_q, const Real weight){
			
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

	}; // class PoissonEvalQuadraturePoint<Geometry, Basis>

} // namespace pdesolver::fem::eval

#endif
