#ifndef POISSON_BOUNDARYVALUEFUNCTION_HPP
#define POISSON_BOUNDARYVALUEFUNCTION_HPP

#include <cmath>
#include <utility>

#include "fem/eval/EvalFunction.hpp"

namespace pdesolver::fem::boundary {

	template<Index SpatialDim, Index numDOFs, class Callable>
	struct PoissonBoundaryValueFunction {

		static constexpr Index NumComponents = numDOFs;

		Callable f;

		constexpr PoissonBoundaryValueFunction(Callable func) : f(std::move(func)) {}

		void eval(const Real time, const Real* x, Real* outValue) const {
			f(time, x, outValue);
		}

		void evalGradient(const Real, const Real*, Real*) const {}

	}; // struct PoissonBoundaryValueFunction

} // namespace pdesolver::fem::boundary

#endif
