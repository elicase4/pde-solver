#ifndef HEATEQUATION_SOURCEFORM_HPP
#define HEATEQUATION_SOURCEFORM_HPP

#include <type_traits>
#include <utility>

#include "fem/form/LinearForm.hpp"
#include "equations/heateq/eval/SourceFunction.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePointVolume, typename SourceFunction>
	requires fem::eval::EvalFunction<SourceFunction>
	struct SourceForm {
		
		SourceFunction sourceFunction;
		
		template<typename... Args>
		requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, SourceForm> && ...)))
		constexpr SourceForm(Args&&... args) : sourceFunction(std::forward<Args>(args)...) {}

		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointVolume& qp, const Real*, Real* Fe) const {
			
			Real val[SourceFunction::NumComponents];
			sourceFunction.eval(qp.time, qp.x, val);

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index i = 0; i < SourceFunction::NumComponents; ++i) {
					Fe[a + i] += (val[i] * qp.N[a]) * qp.measure * qp.w;
				}
			}
		}

	};// struct SourceForm

} // namespace pdesolver::equations::heateq

#endif
