#ifndef PDESOLVER_EVALQUADRATUREPOINT_HPP
#define PDESOLVER_EVALQUADRATUREPOINT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename QuadraturePoint, typename Element>
			concept EvalQuadraturePoint = requires(QuadraturePoint qp, const Element& elem, const Real* xi, const Real w) {

				{ qp.evaluate(elem, xi, w) } -> std::same_as<void>;

			}; // concept EvalElement

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
