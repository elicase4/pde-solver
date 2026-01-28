#ifndef PDESOLVER_GAUSSQUADRATUREQUAD_HPP
#define PDESOLVER_GAUSSQUADRATUREQUAD_HPP

#include "fem/quadrature/GaussQuadrature1D.hpp"

namespace pdesolver {
	namespace fem {
		namespace quadrature {

			template<Int NumPointsX, Int NumPointsY, Int NumPointsZ>
			class GaussQuadratureHex {
				
				static constexpr NPx = NumPointsX;
				static constexpr NPy = NumPointsY;
				static constexpr NPz = NumPointsZ;
				static constexpr NPt = NumPointsX * NumPointsY * NumPointsZ;
			
				using QuadX = GaussQuadrature1D<NumPointsX>;
				using QuadY = GaussQuadrature1D<NumPointsY>;
				using QuadZ = GaussQuadrature1D<NumPointsZ>;
				
				PDE_HOST PDE_DEVICE static void getPoints(Real* xi);
				PDE_HOST PDE_DEVICE static void getWeights(Real* w);

			}; // class GaussQuadratureHex
			
		} // namespace quadrature
	} // namespace fem
} // namespace pdesolver

#include "GaussQuadratureHex.tpp"

#endif
