#ifndef PDESOLVER_EVALMODEL_HPP
#define PDESOLVER_EVALMODEL_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include <concepts>

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Model>
			concept EvalModel = requires (const Model m, const Real* xi, Real* outValue, Real* outGrad) {

				{ m.eval(xi, outValue) };
				{ m.evalGradient(xi, outGrad) };

			}; // concept EvalModel

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
