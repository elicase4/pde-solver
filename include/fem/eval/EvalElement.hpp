#ifndef PDESOLVER_EVALELEMENT_HPP
#define PDESOLVER_EVALELEMENT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Element>
			concept EvalElement = requires(Element elem, const Real* nodeCoords, const Real time) {

				{ Element::SpatialDim } -> std::convertible_to<Int>;
				{ Element::ParametricDim } -> std::convertible_to<Int>;
				{ Element::NumNodes } -> std::convertible_to<Index>;

				{ elem.bind(nodeCoords, time) } -> std::same_as<void>;
				{ elem.quadLoop([](const auto& qp) {}) } -> std::same_as<void>;

			}; // concept EvalElement

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
