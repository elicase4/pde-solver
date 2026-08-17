#ifndef PDESOLVER_EVALQUADRATUREPOINTBOUNDARY_HPP
#define PDESOLVER_EVALQUADRATUREPOINTBOUNDARY_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename QuadraturePointBoundary, typename Element, typename Basis, typename Geometry>
			concept EvalQuadraturePointBoundary = requires(const QuadraturePointBoundary qp, const Element& elem, const Int rngID, const Real* xi, const Real w) {

				{ QuadraturePointBoundary::SpatialDim } -> std::convertible_to<Index>;
				{ QuadraturePointBoundary::ParametricDim } -> std::convertible_to<Index>;

				{ qp.nodesPerElement() } -> std::convertible_to<Index>;
				{ qp.nodesPerFace() } -> std::convertible_to<Index>;
				{ qp.evaluate(xi, w) } -> std::same_as<void>;

				// face-local node data -- always index with a face-local a=0..nodesPerFace()-1
				{ qp.faceNodeLocalIDs[0] } -> std::convertible_to<Index>;
				{ qp.Nface[0] } -> std::convertible_to<Real>;

				{ qp.normal[0] } -> std::convertible_to<Real>;
				{ qp.x[0] } -> std::convertible_to<Real>;
				{ qp.time } -> std::convertible_to<Real>;
				{ qp.w } -> std::convertible_to<Real>;

			}; // concept EvalQuadraturePointBoundary

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
