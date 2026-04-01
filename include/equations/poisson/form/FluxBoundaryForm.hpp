#ifndef POISSON_FLUXBOUNDARYFORM_HPP
#define POISSON_FLUXBOUNDARYFORM_HPP

#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/eval/EvalQuadraturePointBoundary.hpp"
#include "fem/form/LinearForm.hpp"

namespace pdesolver::fem::form {

	template<typename QuadraturePointBoundary, typename FluxFunction>
	struct PoissonFluxBoundaryForm {

		FluxFunction fluxFunction;
		constexpr PoissonFluxBoundaryForm(FluxFunction src) : fluxFunction(std::move(src)) {}
		
		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointBoundary& qp, const Real*, Real* Fe) const {
		
			Real val[FluxFunction::SpatialDim * FluxFunction::NumComponents];
			fluxFunction.eval(qp.time, qp.x, val);

			// element vector assembly contribution
			for (Index i = 0; i < FluxFunction::NumComponents; ++i) {
				for (Index a = 0; a < qp.NodesPerFace; ++a){
					for (Index sD = 0; sD < FluxFunction::SpatialDim; ++sD){
						Fe[a * FluxFunction::NumComponents + i] += (val[i*FluxFunction::SpatialDim + sD] * qp.N[a]) * qp.normal[sD] * qp.w;
					}
				}
			}
		}

	};// struct PoissonFluxBoundaryForm

} // namespace pdesolver::fem::form

#endif
