#ifndef POISSON_SOURCEFORM_HPP
#define POISSON_SOURCEFORM_HPP

#include "fem/form/LinearForm.hpp"
#include "equations/poisson/eval/SourceFunction.hpp"

namespace pdesolver::fem::form {

	template<typename QuadraturePointVolume, typename SourceFunction>
	requires fem::eval::EvalFunction<SourceFunction>
	struct PoissonSourceForm {
		
		SourceFunction sourceFunction;
		constexpr PoissonSourceForm(SourceFunction src) : sourceFunction(std::move(src)) {}

		PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointVolume& qp, const Real*, Real* Fe) const {
			
			Real val[SourceFunction::NumComponents];
			sourceFunction.eval(qp.time, qp.x, val);

			// element vector assembly contribution
			for (Index i = 0; i < SourceFunction::NumComponents; ++i) {
				for (Index a = 0; a < qp.NodesPerElement; ++a){
					Fe[a*SourceFunction::NumComponents + i] += (val[i] * qp.N[a]) * qp.measure * qp.w;
				}
			}
		}

	};// struct PoissonSourceForm

} // namespace pdesolver::fem::form

#endif
