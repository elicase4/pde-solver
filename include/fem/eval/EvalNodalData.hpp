#ifndef PDESOLVER_EVALNODALDATA_HPP
#define PDESOLVER_EVALNODALDATA_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			// node-indexed, no interpolation -- distinct from EvalField, which interpolates an
			// already-gathered element-local buffer at a quadrature point
			template<typename NodalData>
			concept EvalNodalData = requires (const NodalData d, Index nodeID, Real* outValue) {

				{ NodalData::NumComponents } -> std::convertible_to<Index>;
				{ d.eval(nodeID, outValue) } -> std::same_as<void>;

			}; // concept EvalNodalData

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
