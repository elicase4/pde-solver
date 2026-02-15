#ifndef PDESOLVER_EVALELEMENT_HPP
#define PDESOLVER_EVALELEMENT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Element>
			concept EvalElement = requires(Element elem, const Real* nodeCoords, const Real* xi) {
				{ Element::SpatialDim } -> std::convertible_to<Int>;
				{ Element::ParametricDim } -> std::convertible_to<Int>;
				{ Element::NumNodes }  -> std::convertible_to<Int>;

				{ elem.bindElement(nodeCoords) } -> std::same_as<void>;
				{ elem.evaluate(xi) } -> std::same_as<void>;

			}; // concept EvalElement

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
