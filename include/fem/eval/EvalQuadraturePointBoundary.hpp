#ifndef PDESOLVER_EVALQUADRATUREPOINTBOUNDARY_HPP
#define PDESOLVER_EVALQUADRATUREPOINTBOUNDARY_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

#include "fem/boundary/BoundaryCondition.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename QuadraturePointBoundary, typename Element, typename Function>
			concept EvalQuadraturePointBoundary = requires(const QuadraturePointBoundary qp, const Element& elem, const BoundaryCondition<Function>& bc, const Int rngID, const Real* xi, const Real w) {

				{ QuadraturePointBoundary::SpatialDim } -> std::convertible_to<Index>;
				{ QuadraturePointBoundary::ParametricDim } -> std::convertible_to<Index>;
				{ QuadraturePointBoundary::NodesPerFace } -> std::convertible_to<Index>;

				{ qp.evaluate(xi, w) } -> std::same_as<void>;

			}; // concept EvalQuadraturePointBoundary

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
