#ifndef PDESOLVER_EVALELEMENT_HPP
#define PDESOLVER_EVALELEMENT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename Ctx>
			concept EvalElement = requires(Ctx ctx, const Real* nodeCoords, const Real* xi) {
				{ Ctx::SpatialDim } -> std::convertible_to<Int>;
				{ Ctx::ParametricDim } -> std::convertible_to<Int>;
				{ Ctx::NumNodes }  -> std::convertible_to<Int>;

				{ ctx.bindElement(nodeCoords) } -> std::same_as<void>;
				{ ctx.evaluate(xi) } -> std::same_as<void>;
				{ ctx.detJ } -> std::convertible_to<Real>;

			}; // concept EvalElement

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
