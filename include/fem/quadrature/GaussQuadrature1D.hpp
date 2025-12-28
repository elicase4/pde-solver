#ifndef PDESOLVER_GAUSSQUADRATURE1D_HPP
#define PDESOLVER_GAUSSQUADRATURE1D_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace quadrature {
			
			template<Int NumPoints>
			class GaussQuadrature1D {
				PDE_HOST PDE_DEVICE static void getPoints(Real* xi);
				PDE_HOST PDE_DEVICE static void getWeights(Real* w);
			}; // class GaussQuadrature1D
			
			// 1-point implementation
			template <>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<1>::getPoints(Real* xi){
				xi[0] = 0.0;
			}

			template <>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<1>::getWeights(Real* w){
				w[0] = 2.0;
			}
			
			// 2-point implementation
			template<>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<2>::getPoints(Real* xi) {
				xi[0] = -0.5773502691896257;  // -1/sqrt(3)
				xi[1] =  0.5773502691896257;  //  1/sqrt(3)
			}

			template<>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<2>::getWeights(Real* w) {
				w[0] = 1.0;
				w[1] = 1.0;
			}

			// 3-point implementation
			template<>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<3>::getPoints(Real* xi) {
				xi[0] = -0.7745966692414834;  // -sqrt(3/5)
				xi[1] =  0.0;
				xi[2] =  0.7745966692414834;  //  sqrt(3/5)
			}

			template<>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<3>::getWeights(Real* w) {
				w[0] = 0.5555555555555556;  // 5/9
				w[1] = 0.8888888888888889;  // 8/9
				w[2] = 0.5555555555555556;  // 5/9
			}

			// 4-point implementation
			template<>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<4>::getPoints(Real* xi) {
				xi[0] = -0.8611363115940526;
				xi[1] = -0.3399810435848563;
				xi[2] =  0.3399810435848563;
				xi[3] =  0.8611363115940526;
			}

			template<>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<4>::getWeights(Real* w) {
				w[0] = 0.3478548451374538;
				w[1] = 0.6521451548625461;
				w[2] = 0.6521451548625461;
				w[3] = 0.3478548451374538;
			}

			// 5-point implementation
			template<>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<5>::getPoints(Real* xi) {
				xi[0] = -0.9061798459386640;
				xi[1] = -0.5384693101056831;
				xi[2] =  0.0;
				xi[3] =  0.5384693101056831;
				xi[4] =  0.9061798459386640;
			}

			template<>
			PDE_HOST PDE_DEVICE void GaussQuadrature1D<5>::getWeights(Real* w) {
				w[0] = 0.2369268850561891;
				w[1] = 0.4786286704993665;
				w[2] = 0.5688888888888889;
				w[3] = 0.4786286704993665;
				w[4] = 0.2369268850561891;
			}

		} // namespace quadrature
	} // namespace fem
} // namespace pdesolver

#endif
