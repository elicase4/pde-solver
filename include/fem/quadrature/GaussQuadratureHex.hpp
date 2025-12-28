#ifndef PDESOLVER_GAUSSQUADRATUREQUAD_HPP
#define PDESOLVER_GAUSSQUADRATUREQUAD_HPP

#include "fem/quadrature/GaussQuadrature1D.hpp"

namespace pdesolver {
	namespace fem {
		namespace quadrature {

			template<Int NumPointsX, Int NumPointsY, Int NumPointsZ>
			class GaussQuadratureHex {
				
				using QuadX = GaussQuadrature1D<NumPointsX>;
				using QuadY = GaussQuadrature1D<NumPointsY>;
				using QuadZ = GaussQuadrature1D<NumPointsZ>;
				
				PDE_HOST PDE_DEVICE static void getPoints(Real* xi);
				PDE_HOST PDE_DEVICE static void getWeights(Real* w);
			
			}; // class GaussQuadratureHex
			
			// implementation points
			template<Int NumPointsX, Int NumPointsY, Int NumPointsZ>
			PDE_HOST PDE_DEVICE void GaussQuadratureHex<NumPointsX, NumPointsY, NumPointsZ>::getPoints(Real* xi){
				
				Real xi_xi[NumPointsX];
				Real xi_eta[NumPointsY];
				Real xi_zeta[NumPointsZ];

				QuadX::getPoints(xi_xi);
				QuadY::getPoints(xi_eta);
				QuadZ::getPoints(xi_zeta);

				// compute tensor product
				Index q;
				for (Index k = 0; k < NumPointsZ; ++k) {
					for (Index j = 0; j < NumPointsY; ++j) {
						for (Index i = 0; i < NumPointsX; ++i) {
							xi[3*q]   = xi_xi[i];
							xi[3*q+1] = xi_eta[j];
							xi[3*q+2] = xi_zeta[k];
							q++;
						}
					}
				}
			}

			// implementation weights
			template<Int NumPointsX, Int NumPointsY, Int NumPointsZ>
			PDE_HOST PDE_DEVICE void GaussQuadratureQuad<NumPointsX, NumPointsY, NumPointsZ>::getWeights(Real* w){
				
				Real w_xi[NumPointsX];
				Real w_eta[NumPointsY];
				Real w_zeta[NumPointsZ];

				QuadX::getPoints(w_xi);
				QuadY::getPoints(w_eta);
				QuadZ::getPoints(w_zeta);

				// compute tensor product
				Index q;
				for (Index k = 0; k < NumPointsZ; ++k) {
					for (Index j = 0; j < NumPointsY; ++j) {
						for (Index i = 0; i < NumPointsX; ++i) {
							w[q] = w_xi[i] * w_eta[j] * w_zeta[k];
							q++;
						}
					}
				}
			}

		} // namespace quadrature
	} // namespace fem
} // namespace pdesolver

#endif
