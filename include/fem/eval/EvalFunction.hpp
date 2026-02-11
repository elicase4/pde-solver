#ifndef PDESOLVER_EVALFUNCTION_HPP
#define PDESOLVER_EVALFUNCTION_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Func>
			concept EvalFunction = requires (const Func f, const Real t, const Real* x, Real* out) {
				
				{ Func::NumComponents() } -> std::convertible_to<Index>;
				
				{ f.evaluate(t, x, out) } -> std::same_as<void>;
			
			}; // concept EvalFunction

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
