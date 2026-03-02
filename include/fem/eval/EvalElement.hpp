#ifndef PDESOLVER_EVALELEMENT_HPP
#define PDESOLVER_EVALELEMENT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Element>
			concept EvalElement = requires(Element elem, const Real* nodeCoords, const Real time) {

				{ Element::SpatialDim } -> std::convertible_to<Index>;
				{ Element::ParametricDim } -> std::convertible_to<Index>;
				{ Element::NodesPerElement } -> std::convertible_to<Index>;
				
				{ elem.bindElement(nodeCoords, time) } -> std::same_as<void>;

			}; // concept EvalElement

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
