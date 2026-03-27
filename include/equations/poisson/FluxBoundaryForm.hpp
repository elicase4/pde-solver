#ifndef POISSON_FLUXBOUNDARYFORM_HPP
#define POISSON_FLUXBOUNDARYFORM_HPP

#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/eval/EvalQuadraturePointBoundary.hpp"
#include "fem/form/LinearForm.hpp"

namespace pdesolver::fem::form {

	template<typename QuadraturePointBoundary, typename BoundaryCondition, typename Function, Int SpatialDim>
	struct PoissonFluxBoundaryForm {
		
		Function source;

		constexpr PoissonFluxBoundaryForm(BoundaryCondition bc) : source(std::move(bc.f)) {}

		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointBoundary& qp, const Real*, Real* Fe) const {
			
			Real val[1];
			source.eval(qp.time, qp.x, val);

			// element vector assembly contribution
			for (Index a = 0; a < qp.NodesPerFace; ++a){
				Fe[a] += (val[0] * qp.N[a]) * qp.normal * qp.w;
			}
		}

	};// struct PoissonFluxBoundaryForm

} // namespace pdesolver::fem::form

#endif
