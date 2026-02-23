#ifndef PDESOLVER_GAUSSQUADRATURE1D_HPP
#define PDESOLVER_GAUSSQUADRATURE1D_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace quadrature {
			
			template<Int NumPoints>
			class GaussQuadrature1D {
				
				static constexpr Int NP = NumPoints;

				PDE_HOST PDE_DEVICE static void getPoints(Real* xi);
				PDE_HOST PDE_DEVICE static void getWeights(Real* w);

			}; // class GaussQuadrature1D
			
		} // namespace quadrature
	} // namespace fem
} // namespace pdesolver

#include "GaussQuadrature1D.tpp"

#endif
