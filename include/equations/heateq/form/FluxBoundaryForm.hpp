#ifndef HEATEQUATION_FLUXBOUNDARYFORM_HPP
#define HEATEQUATION_FLUXBOUNDARYFORM_HPP

#include <type_traits>
#include <utility>

#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/eval/EvalQuadraturePointBoundary.hpp"
#include "fem/form/LinearForm.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePointBoundary, typename FluxFunction>
	struct FluxBoundaryForm {

		FluxFunction fluxFunction;
		
		template<typename... Args>
		requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, FluxBoundaryForm> && ...)))
		constexpr FluxBoundaryForm(Args&&... args) : fluxFunction(std::forward<Args>(args)...) {}
		
		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointBoundary& qp, const Real*, Real* Fe) const {
		
			Real val[FluxFunction::SpatialDim * FluxFunction::NumComponents];
			fluxFunction.eval(qp.time, qp.x, val);

			for (Index a = 0; a < qp.nodesPerFace(); ++a){
				for (Index i = 0; i < FluxFunction::NumComponents; ++i) {
					for (Index sD = 0; sD < FluxFunction::SpatialDim; ++sD){
						Fe[a*FluxFunction::NumComponents + i] += (val[i*FluxFunction::SpatialDim + sD] * qp.Nface[a]) * qp.normal[sD] * qp.w;
					}
				}
			}

		}

	};// struct FluxBoundaryForm

} // namespace pdesolver::equations::heateq

#endif
