#ifndef PDESOLVER_EVALFIELD_HPP
#define PDESOLVER_EVALFIELD_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Field>
			concept EvalField = requires (const Real t, const Real* x, Real* v) {
				{ EvalField::evaluate(t, x, v) } -> std::same_as<void>;
				{ EvalField::numComponents() } -> std::convertible_to<Index>;
			}; // concept EvalField

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
