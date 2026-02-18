#ifndef POISSON_EVALQUADRATUREPOINT_HPP
#define POISSON_EVALQUADRATUREPOINT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/EvalGeometry.hpp"
#include "fem/eval/EvalQuadraturePoint.hpp"
#include "fem/eval/EvalModel.hpp"

namespace pdesolver::fem::eval {

	template<typename Geometry, typename Basis, typename ConductivityModel>
	class PoissonEvalQuadraturePoint<Geometry, Basis, Geometry, ConductivityModel> {
	public:
		
		// physical coordinate
		Real x[Geometry::SpatialDim];

		// quadrature
		Real xi[Geometry::ParametricDim];
		Real w;

		// ref basis values
		Real N[Geometry::NumNodes];

		// ref gradients
		Real dNdxi[Geometry::ParametricDim*Geometry::NumNodes];

		// physical gradients
		Real dNdx[Geometry::SpatialDim*Geometry::NumNodes];

		// geometry
		Real J[Geometry::SpatialDim*Geometry::ParametricDim];
		Real g[Geometry::ParametricDim*Geometry::ParametricDim];

		// measure
		Real measure;

		// conductivity coefficient
		Real K[Geometry::SpatialDim * Geometry::SpatialDim];

		// rhs function
		Real rhsF[Geometry::SpatialDim];

		PDE_HOST PDE_DEVICE evaluate(const Real* coords, const Real* xi_q, const Real weight){
			
			// set quad info
			for (Index pD = 0; pD < ParametricDim; ++pD){
				xi[pD] = xi_q[pD];
			}
			w = weight;
			
			// evaluate basis
			Basis::evaluate(xi, N);
			Basis::evaluateGradient(xi, dNdxi);

			// geometry
			Geometry::mapToPhysical(coords, N, x);
			Geoemtry::computeJacobian(coords, dNdxi, J);
			Geoemtry::computeMetric(J, g);
			measure = Geometry::computeMeasure(g);

			// transforms
			Geometry::transformGradient(J, g, dNdxi, dNdx);

			// evalaute conductivity model
			ConductivityModel::eval(qp);
			
		}

	}; // class PoissonEvalQuadraturePoint<Geometry, Basis, ConductivityModel, SourceTerm>

} // namespace pdesolver::fem::eval

#endif
