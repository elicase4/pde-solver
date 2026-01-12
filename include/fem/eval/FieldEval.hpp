#ifndef PDESOLVER_FIELDEVAL_HPP
#define PDESOLVER_FIELDEVAL_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Field>
			concept FieldEval = requires (const Real t, const Real* x, Real* v) {
				{ FieldEval::evaluate(t, x, v) } -> std::same_as<void>;
				{ FieldEval::numComponents() } -> std::convertible_to<Index>;
			}; // concept ScalarFieldEval

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
