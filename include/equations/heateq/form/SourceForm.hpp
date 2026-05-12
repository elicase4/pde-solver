#ifndef POISSON_SOURCEFORM_HPP
#define POISSON_SOURCEFORM_HPP

#include "fem/form/LinearForm.hpp"
#include "equations/heateq/eval/SourceFunction.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePointVolume, typename SourceFunction>
	requires fem::eval::EvalFunction<SourceFunction>
	struct SourceForm {
		
		SourceFunction sourceFunction;
		constexpr SourceForm(SourceFunction src) : sourceFunction(std::move(src)) {}

		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointVolume& qp, const Real*, Real* Fe) const {
			
			Real val[SourceFunction::NumComponents];
			sourceFunction.eval(qp.time, qp.x, val);

			// element vector assembly contribution
			for (Index a = 0; a < qp.NodesPerElement; ++a){
				for (Index i = 0; i < SourceFunction::NumComponents; ++i) {
					Fe[a + i] += (val[i] * qp.N[a]) * qp.measure * qp.w;
				}
			}
		}

	};// struct SourceForm

} // namespace pdesolver::equations::heateq

#endif
