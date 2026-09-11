#ifndef HEATEQUATION_QUANTITY_HEATFLUXINTEGRAND_HPP
#define HEATEQUATION_QUANTITY_HEATFLUXINTEGRAND_HPP

#include "config/Platform.hpp"
#include "core/Types.hpp"

#include "equations/heateq/eval/EvalField.hpp"

namespace pdesolver::equations::heateq::quantity {

	template<typename QuadraturePointBoundary>
	class HeatFluxIntegrand {
	public:

		static constexpr Index NumComponents = 1;

		static constexpr const char* BaseUnit = "W";

		PDE_HOST PDE_DEVICE void computeElementLevelValue(const QuadraturePointBoundary& qp, const Real* Ue, Real* out) const {

			static constexpr Index SpatialDim = QuadraturePointBoundary::SpatialDim;

			Real gradT[SpatialDim];
			EvalField().evalGradient(qp, Ue, gradT);

			// flux = -(K*gradT) . normal
			Real flux = 0;
			for (Index i = 0; i < SpatialDim; ++i) {
				Real KgradT_i = 0;
				for (Index j = 0; j < SpatialDim; ++j) {
					KgradT_i += qp.K[i*SpatialDim + j] * gradT[j];
				}
				flux -= KgradT_i * qp.normal[i];
			}

			out[0] = flux * qp.w;

		}

	}; // class HeatFluxIntegrand

} // namespace pdesolver::equations::heateq::quantity

#endif
