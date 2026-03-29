#ifndef POISSON_BOUNDARYFLUXFUNCTION_HPP
#define POISSON_BOUNDARYFLUXFUNCTION_HPP

#include <cmath>
#include <utility>

#include "core/Types.hpp"
#include "fem/eval/EvalFunction.hpp"

namespace pdesolver::fem::boundary {

	template<Int SpatialDim, class Callable>
	struct PoissonBoundaryFluxFunction {

		static constexpr Index NumComponents = SpatialDim;

		Callable f;

		constexpr PoissonBoundaryFluxFunction(Callable func) : f(std::move(func)) {}

		void eval(const Real time, const Real* x, Real* outValue) const {
			f(time, x, outValue);
		}

		void evalGradient(const Real, const Real*, Real*) const {}

	}; // struct PoissonBoundaryFluxFunction

} // namespace pdesolver::fem::boundary

#endif
