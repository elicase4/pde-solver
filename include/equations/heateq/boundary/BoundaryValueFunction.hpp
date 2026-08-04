#ifndef HEATEQUATION_BOUNDARYVALUEFUNCTION_HPP
#define HEATEQUATION_BOUNDARYVALUEFUNCTION_HPP

#include <cmath>
#include <type_traits>
#include <utility>

#include "fem/eval/EvalFunction.hpp"

namespace pdesolver::equations::heateq {

	template<Index SpatialDimension, Index numDOFs, class Callable>
	struct BoundaryValueFunction {

		static constexpr Index NumComponents = numDOFs;
		static constexpr Index SpatialDim = SpatialDimension;

		Callable f;
		
		template<typename... Args>
		requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, BoundaryValueFunction> && ...)))
		constexpr BoundaryValueFunction(Args&&... args) : f(std::forward<Args>(args)...) {}

		void eval(const Real time, const Real* x, Real* outValue) const {
			f(time, x, outValue);
		}

		void evalGradient(const Real, const Real*, Real*) const {}

	}; // struct BoundaryValueFunction

} // namespace pdesolver::equations::heateq

#endif
