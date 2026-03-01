#ifndef PDESOLVER_EVALELEMENT_HPP
#define PDESOLVER_EVALELEMENT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Element, Index SpatialDim>
			concept EvalElement = requires(Element elem, const Real* nodeCoords, const Real time) {

				{ SpatialDim } -> std::convertible_to<Index>;
				{ Element::NodesPerElement } -> std::convertible_to<Index>;
				
				{ elem.bind(nodeCoords, time) } -> std::same_as<void>;

			}; // concept EvalElement

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
