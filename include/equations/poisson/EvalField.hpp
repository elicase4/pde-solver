#ifndef POISSON_EVALFIELD_HPP
#define POISSON_EVALFIELD_HPP

#include "fem/eval/EvalField.hpp"

namespace pdesolver::fem::eval {

	template<Index NumNodes, Int SpatialDim>
	struct PoissonField {

		static constexpr Index NumComponents = 1;

		PDE_HOST PDE_DEVICE void eval(const quto& qp, const Real* Ue, Real* outValue) const {
		}

		PDE_HOST PDE_DEVICE void evalGradient(const quto& qp, const Real* Ue, Real* outValue) const {
		}

	}; // struct PoissonField

}

#endif
