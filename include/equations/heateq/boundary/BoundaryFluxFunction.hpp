#ifndef POISSON_BOUNDARYFLUXFUNCTION_HPP
#define POISSON_BOUNDARYFLUXFUNCTION_HPP

#include <cmath>
#include <utility>

#include "core/Types.hpp"
#include "fem/eval/EvalFunction.hpp"

namespace pdesolver::equations::heateq {

	template<Index SpatialDimension, Index numDOFs, class Callable>
	struct BoundaryFluxFunction {

		static constexpr Index NumComponents = numDOFs;
		static constexpr Index SpatialDim = SpatialDimension;

		Callable f;

		constexpr BoundaryFluxFunction(Callable func) : f(std::move(func)) {}

		void eval(const Real time, const Real* x, Real* outValue) const {
			f(time, x, outValue);
		}

		void evalGradient(const Real, const Real*, Real*) const {}

	}; // struct BoundaryFluxFunction

} // namespace pdesolver::equations::heateq

#endif
