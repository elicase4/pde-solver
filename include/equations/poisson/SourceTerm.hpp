#ifndef POISSON_SOURCETERM_HPP
#define POISSON_SOURCETERM_HPP

#include <cmath>
#include <utility>

#include "fem/eval/EvalFunction.hpp"

namespace pdesolver::fem::eval {

	template<Int SpatialDim, class Callable>
	struct PoissonSourceTerm {

		static constexpr Index NumComponents = 1;

		Callable f;

		constexpr PoissonSourceTerm(Callable func) : f(std::move(func)) {}

		void eval(const Real time, const Real* x, Real* outValue) const {
			outValue[0] = f(time, x);
		}

	}; // struct PoissonSourceTerm

} // namespace pdesolver::fem::eval

#endif
