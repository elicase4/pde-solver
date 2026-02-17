#ifndef PDESOLVER_EVALQUADRATUREPOINT_HPP
#define PDESOLVER_EVALQUADRATUREPOINT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename QuadraturePoint>
			concept EvalQuadraturePoint = requires(QuadraturePoint qp, const Real* NodeCoords, const Real* xi, const Real w) {

				{ qp.evaluate(nodeCoords, xi, w) } -> std::same_as<void>;

			}; // concept EvalQuadraturePoint

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
