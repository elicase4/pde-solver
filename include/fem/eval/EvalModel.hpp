#ifndef PDESOLVER_EVALMODEL_HPP
#define PDESOLVER_EVALMODEL_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename M>
			concept EvalModel = requires (const M m, const typename M::Context& ctx) {

				{ m.evaluate(ctx) };
				{ m.derivative(ctx) };

			}; // concept EvalModel

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
