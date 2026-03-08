#ifndef POISSON_LINEARFORM_HPP
#define POISSON_LINEARFORM_HPP

#include "fem/form/LinearForm.hpp"
#include "equations/poisson/SourceTerm.hpp"

namespace pdesolver::fem::form {

	template<typename QuadraturePoint, Int SpatialDim, typename SourceFunction>
	requires fem::eval::EvalFunction<SourceFunction>
	struct PoissonLinearForm {
		
		SourceFunction source;
		constexpr PoissonLinearForm(SourceFunction src) : source(std::move(src)) {}

		PDE_HOST PDE_DEVICE void computeElementLevel(const QuadraturePoint& qp, const Real*, Real* Fe) const {
			
			Real val[1];
			source.eval(qp.time, qp.x, val);

			// element vector assembly contribution
			for (Index a = 0; a < qp.NodesPerElement; ++a){
				Fe[a] += (val[0] * qp.N[a]) * qp.measure * qp.w;
			}
		}

	};// struct PoissonLinearForm

} // namespace pdesolver::fem::form

#endif
