#ifndef PDESOLVER_GAUSSQUADRATUREQUAD_HPP
#define PDESOLVER_GAUSSQUADRATUREQUAD_HPP

#include "fem/quadrature/GaussQuadrature1D.hpp"

namespace pdesolver {
	namespace fem {
		namespace quadrature {

			template<Index NumPointsX, Index NumPointsY, Index NumPointsZ>
			class GaussQuadratureHex {
			public:

				static constexpr Index NumPointsXi = NumPointsX;
				static constexpr Index NumPointsEta = NumPointsY;
				static constexpr Index NumPointsZeta = NumPointsZ;
				static constexpr Index NumPointsTotal = NumPointsX * NumPointsY * NumPointsZ;
			
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
