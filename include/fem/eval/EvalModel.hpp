#ifndef PDESOLVER_EVALFIELD_HPP
#define PDESOLVER_EVALFIELD_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>


namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename M>
			concept EvalModel = requires (const M m, const typename M::Context& ctx) {

				{ m.evaluate(ctx) } -> std::same_as<void>;
				{ m.derivatives(ctx) } -> std::same_as<void>;

			}; // concept EvalModel

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
