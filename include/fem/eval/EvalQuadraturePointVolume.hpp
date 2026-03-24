#ifndef PDESOLVER_EVALQUADRATUREPOINTVOLUME_HPP
#define PDESOLVER_EVALQUADRATUREPOINTVOLUME_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename QuadraturePointVolume, typename Element>
			concept EvalQuadraturePointVolume = requires(const QuadraturePointVolume qp, const Element& elem, const Real* xi, const Real w) {

				{ QuadraturePointVolume::SpatialDim } -> std::convertible_to<Index>;
				{ QuadraturePointVolume::ParametricDim } -> std::convertible_to<Index>;
				{ QuadraturePointVolume::NodesPerElement } -> std::convertible_to<Index>;

				{ qp.evaluate(xi, w) } -> std::same_as<void>;

			}; // concept EvalQuadraturePointVolume

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
