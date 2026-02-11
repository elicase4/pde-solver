#ifndef PDESOLVER_EVALFIELD_HPP
#define PDESOLVER_EVALFIELD_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename F>
			concept EvalField = requires (const F f, Int qp) {

				{ static constexpr F::NumComponents() } -> std::convertible_to<Index>;

				{ f.value(qp) } -> std::same_as<const Real*>;
				{ f.grad(qp) } -> std::same_as<const Real*>;

			}; // concept EvalField

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
