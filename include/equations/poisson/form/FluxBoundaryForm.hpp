#ifndef POISSON_FLUXBOUNDARYFORM_HPP
#define POISSON_FLUXBOUNDARYFORM_HPP

#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/eval/EvalQuadraturePointBoundary.hpp"
#include "fem/form/LinearForm.hpp"

namespace pdesolver::fem::form {

	template<typename QuadraturePointBoundary, Int SpatialDim>
	struct PoissonFluxBoundaryForm {
		
		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointBoundary& qp, const Real*, Real* Fe) const {
			
			// element vector assembly contribution
			for (Index a = 0; a < qp.NodesPerFace; ++a){
				for (Index i = 0; i < SpatialDim; ++i){
					Fe[a] += (qp.fluxVal[i] * qp.N[a]) * qp.normal[i] * qp.w;
				}
			}
		}

	};// struct PoissonFluxBoundaryForm

} // namespace pdesolver::fem::form

#endif
