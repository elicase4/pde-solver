#ifndef PDESOLVER_EVALMODEL_HPP
#define PDESOLVER_EVALMODEL_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Model, typename QuadraturePoint>
			concept EvalModel = requires (const Model m, QuadraturePoint& qp) {

				{ m.eval(qp) } -> std::same_as<void>;
				{ m.evalGradient(qp) } -> std::same_as<void>;

			}; // concept EvalModel

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
