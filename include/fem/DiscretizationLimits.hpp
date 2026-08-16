#ifndef PDESOLVER_FEM_DISCRETIZATIONLIMITS_HPP
#define PDESOLVER_FEM_DISCRETIZATIONLIMITS_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace fem {

		// Upper bounds on runtime discretization parameters (basis order, quadrature
		// point count), used to size stack buffers that would otherwise need a
		// template parameter per order/point-count. Raise these if fem::basis /
		// fem::quadrature grow support for a higher order or a finer rule -- they
		// must stay in sync with what Lagrange1D / GaussQuadrature1D actually
		// implement.
		inline constexpr Index kMaxBasisOrder = 3;
		inline constexpr Index kMaxQuadraturePoints1D = 5;

		namespace detail {

			constexpr Index ipow(Index base, Index exp) {
				Index result = 1;
				for (Index i = 0; i < exp; ++i) {
					result *= base;
				}
				return result;
			}

		} // namespace detail

		// Max nodes in a tensor-product element of parametric dimension NPD, at
		// kMaxBasisOrder. NPD stays a compile-time template parameter throughout the
		// codebase (see mesh::ElementFamily / fem::dispatch), so callers already have
		// it in scope wherever a buffer needs sizing.
		template<Index NPD>
		inline constexpr Index kMaxNodesPerElement = detail::ipow(kMaxBasisOrder + 1, NPD);

		// Max total quadrature points in a tensor-product NPD-dimensional volume rule.
		template<Index NPD>
		inline constexpr Index kMaxQuadraturePointsTotal = detail::ipow(kMaxQuadraturePoints1D, NPD);

		// Max total quadrature points on an (NPD-1)-dimensional boundary/face rule.
		template<Index NPD>
		inline constexpr Index kMaxQuadraturePointsTotalBoundary = detail::ipow(kMaxQuadraturePoints1D, NPD - 1);

	} // namespace fem
} // namespace pdesolver

#endif
