#ifndef POISSON_EVALELEMENT_HPP
#define POISSON_EVALELEMENT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/EvalElement.hpp"
#include "fem/eval/EvalModel.hpp"

namespace pdesolver::fem::eval {

	template<typename Geometry, typename Basis>
	class PoissonEvalElement<Geometry, Basis, ConductivityModel> {
	public:
		static constexpr Int ParametricDim = Basis::ParametricDim;
		static constexpr Int SpatialDim = Basis::SpatialDim;
		static constexpr Int NumNodes = Basis::NumNodes;
		
		// node coordinates
		const Real nodeCoords[SpatialDim*NumNodes];
		
		// physical coordinate
		Real x[SpatialDim*NumNodes];

		// time coordinate
		Real t;

		// quadrature
		Real xi[ParametricDim];
		Real w;

		// geometry
		Real J[SpatialDim*ParametricDim];
		Real g[ParametricDim*ParametricDim];

		// basis values
		Real N[NumNodes];

		// ref gradients
		Real dNdxi[ParametricDim*NumNodes];

		// physical gradients
		Real dNdx[SpatialDim*NumNodes];

		// measure
		Real measure;

		// coefficients
		Real K[SpatialDim * SpatialDim];
		Real dK[SpatialDim * SpatialDim];

		PDE_HOST PDE_DEVICE bindElement(const Real* coords, const Real time){
			nodeCoords = coords;
			t = time;
		}

		PDE_HOST PDE_DEVICE evaluate(const Real* xi_q, const Real weight){
			
			// set quad info
			for (Index pD = 0; pD < ParametricDim; ++pD){
				xi[pD] = xi_q[pD];
			}
			w = weight;
			
			// evaluate basis
			Basis::evaluate(xi, N);
			Basis::evaluateGradient(xi, dNdxi);

			// geometry
			Geometry::mapToPhysical(nodeCoords, N, x);
			Geoemtry::computeJacobian(nodeCoords, dNdxi, J);
			Geoemtry::computeMetric(J, g);
			measure = Geometry::computeMeasure(g);

			// transforms
			Geometry::transformGradient(J, g, dNdxi, dNdx);
		}

	}; // class PoissonEvalElement<Geometry, Basis>

} // namespace pdesolver::fem::eval

#endif
