#ifndef POISSON_FLUXBOUNDARYFORM_HPP
#define POISSON_FLUXBOUNDARYFORM_HPP

#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/eval/EvalQuadraturePointBoundary.hpp"
#include "fem/form/LinearForm.hpp"

namespace pdesolver::fem::form {

	template<typename QuadraturePointBoundary, typename BoundaryCondition, typename Function, Int SpatialDim>
	struct PoissonFluxBoundaryForm {
		
		Function flux;

		constexpr PoissonFluxBoundaryForm(BoundaryCondition bc) : flux(std::move(bc.f)) {}

		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointBoundary& qp, const Real*, Real* Fe) const {
			
			Real val[SpatialDim];
			flux.eval(qp.time, qp.x, val);

			// element vector assembly contribution
			for (Index a = 0; a < qp.NodesPerFace; ++a){
				for (Index i = 0; i < SpatialDim; ++i){
					Fe[a] += (val[i] * qp.N[a]) * qp.normal[i] * qp.w;
				}
			}
		}

	};// struct PoissonFluxBoundaryForm

} // namespace pdesolver::fem::form

#endif
