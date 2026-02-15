#ifndef PDESOLVER_EVALFUNCTION_HPP
#define PDESOLVER_EVALFUNCTION_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Function>
			concept EvalFunction = requires (const Function f, const Real t, const Real* x, Real* outValue, Real* outGrad) {
				
				{ Function::NumComponents } -> std::convertible_to<Index>;
				
				{ f.eval(t, x, outValue) } -> std::same_as<void>;
				{ f.evalGradient(t, x, outGrad) } -> std::same_as<void>;
			
			}; // concept EvalFunction

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
