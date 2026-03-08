#ifndef PDESOLVER_EVALQUADRATUREPOINT_HPP
#define PDESOLVER_EVALQUADRATUREPOINT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename QuadraturePoint, typename Element>
			concept EvalQuadraturePoint = requires(const QuadraturePoint qp, const Element& elem, const Real* xi, const Real w) {

				{ QuadraturePoint::SpatialDim } -> std::convertible_to<Index>;
				{ QuadraturePoint::ParametricDim } -> std::convertible_to<Index>;
				{ QuadraturePoint::NodesPerElement } -> std::convertible_to<Index>;

				{ qp.evaluate(xi, w) } -> std::same_as<void>;

			}; // concept EvalQuadraturePoint

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
