#ifndef PDESOLVER_EVALFIELD_HPP
#define PDESOLVER_EVALFIELD_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Field>
			concept EvalField = requires (const Field field, const Real* xi, Real* outValue, Real* outGrad) {

				{ Field::NumComponents } -> std::convertible_to<Index>;

				{ field.eval(xi, outValue) } -> std::same_as<void>;
				{ field.evalGradient(xi, outGrad) } -> std::same_as<void>;

			}; // concept EvalField

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
