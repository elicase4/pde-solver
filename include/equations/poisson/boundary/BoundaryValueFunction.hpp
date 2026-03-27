#ifndef POISSON_BOUNDARYVALUEFUNCTION_HPP
#define POISSON_BOUNDARYVALUEFUNCTION_HPP

#include <cmath>
#include <utility>

#include "fem/eval/EvalFunction.hpp"

namespace pdesolver::fem::boundary {

	template<Int SpatialDim, class Callable>
	struct PoissonBoundaryValueFunction {

		static constexpr Index NumComponents = 1;

		Callable f;

		constexpr PoissonBoundaryValueFunction(Callable func) : f(std::move(func)) {}

		void eval(const Real time, const Real* x, Real* outValue) const {
			outValue[0] = f(time, x);
		}

		void evalGradient(const Real, const Real*, Real*) const {}

	}; // struct PoissonBoundaryValueFunction

} // namespace pdesolver::fem::boundary

#endif
