#ifndef HEATEQUATION_EVALFIELD_HPP
#define HEATEQUATION_EVALFIELD_HPP

#include <type_traits>

#include "fem/eval/EvalField.hpp"

namespace pdesolver::equations::heateq {

	// no template params -- node count comes from qp.nodesPerElement() (runtime, matching every
	// other class since the dispatch refactor) and SpatialDim from decltype(qp), so one instance
	// works for both EvalQPVol and EvalQPBdy.
	struct EvalField {

		static constexpr Index NumComponents = 1;

		PDE_HOST PDE_DEVICE void eval(const auto& qp, const Real* Ue, Real* outValue) const {

			Real u = 0;

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				u += qp.N[a] * Ue[a];
			}

			outValue[0] = u;
		}

		PDE_HOST PDE_DEVICE void evalGradient(const auto& qp, const Real* Ue, Real* outGrad) const {

			using QP = std::decay_t<decltype(qp)>;

			for (Index i = 0; i < QP::SpatialDim; ++i){
				outGrad[i] = 0;
			}

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index i = 0; i < QP::SpatialDim; ++i){
					outGrad[i] += qp.dNdx[a*QP::SpatialDim + i] * Ue[a];
				}
			}

		}

	}; // struct EvalField

} // namespace pdesolver::equations::heateq

#endif
