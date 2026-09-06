#ifndef HEATEQUATION_NODALFLUXFORM_HPP
#define HEATEQUATION_NODALFLUXFORM_HPP

#include <utility>

#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/eval/EvalNodalData.hpp"
#include "fem/eval/EvalQuadraturePointBoundary.hpp"
#include "fem/form/LinearForm.hpp"

namespace pdesolver::equations::heateq {

	// data-driven counterpart to FluxBoundaryForm -- instead of evaluating a continuous
	// expression at qp.x, interpolates a per-face-node array (gathered once per face via
	// gatherElementData) at the QP using the same face basis as the test space (matching
	// Galerkin: same basis for data interpolation and test-function weighting).
	template<typename QuadraturePointBoundary, fem::eval::EvalNodalData Source, Index NumComponents_, Index SpatialDim_>
	class NodalFluxForm {
	public:

		static constexpr Index NumComponents = NumComponents_;
		static constexpr Index SpatialDim = SpatialDim_;

		template<typename... Args>
		requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, NodalFluxForm> && ...)))
		constexpr NodalFluxForm(Args&&... args) : source_(std::forward<Args>(args)...) {}

		// per-face, before the quadrature-point loop -- see fem/form/GathersElementData.hpp
		void gatherElementData(const Index* faceNodeGlobalIDs, Index nodesPerFace) const {
			for (Index a = 0; a < nodesPerFace; ++a) {
				source_.eval(faceNodeGlobalIDs[a], &faceDe_[a * NumComponents * SpatialDim]);
			}
		}

		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointBoundary& qp, const Real*, Real* Fe) const {

			Real val[SpatialDim * NumComponents] = {0};

			// interpolate the gathered per-node flux data to this QP
			for (Index b = 0; b < qp.nodesPerFace(); ++b) {
				for (Index i = 0; i < NumComponents; ++i) {
					for (Index sD = 0; sD < SpatialDim; ++sD) {
						val[i*SpatialDim + sD] += qp.Nface[b] * faceDe_[b*NumComponents*SpatialDim + i*SpatialDim + sD];
					}
				}
			}

			// same weak-form accumulation as FluxBoundaryForm
			for (Index a = 0; a < qp.nodesPerFace(); ++a){
				for (Index i = 0; i < NumComponents; ++i) {
					for (Index sD = 0; sD < SpatialDim; ++sD){
						Fe[a*NumComponents + i] += (val[i*SpatialDim + sD] * qp.Nface[a]) * qp.normal[sD] * qp.w;
					}
				}
			}

		}

	private:

		Source source_;
		mutable Real faceDe_[fem::dispatch::kMaxNodesPerElementBoundary<QuadraturePointBoundary::ParametricDim> * NumComponents_ * SpatialDim_];

	}; // class NodalFluxForm

} // namespace pdesolver::equations::heateq

#endif
