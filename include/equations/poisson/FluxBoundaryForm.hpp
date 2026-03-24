#ifndef POISSON_FLUXBOUNDARYFORM_HPP
#define POISSON_FLUXBOUNDARYFORM_HPP

#include "fem/form/LinearForm.hpp"
#include "fem/boundary/BoundaryCondition.hpp"

namespace pdesolver::fem::form {

	template<typename QuadraturePointBoundary, Int SpatialDim, typename BoundaryFunction>
	requires fem::boundary::BoundaryFunction<BFunction>
	struct PoissonFluxBoundaryForm {
		
		BFunction source;
		constexpr PoissonFluxBoundaryForm(BFunction src) : source(std::move(src)) {}

		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointVolume& qp, const Real*, Real* Fe) const {
			
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
