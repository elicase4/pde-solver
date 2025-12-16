#ifndef PDESOLVER_LAGRANGEHEX_HPP
#define PDESOLVER_LAGRANGEHEX_HPP

#include "fem/Lagrange1D.hpp"

namespace pdesolver {
	
	namespace fem {
		
		template <int Px, int Py, int Pz>
		class LagrangeQuad {
		public:
			static constexpr int nodesPerElement = (Px + 1) * (Py + 1) * (Pz + 1);

			static void eval(const Real* xi, Real* N);

			static void evalGradient(const Real* xi, Real* dNdxi, Real* dNdtheta, Real* dNdzeta);

			static void evalHessian(const Real* xi, Real* d2Nd2xi, Real* d2Nd2theta, Real* d2Nd2zeta, 
													Real* d2Ndthetadxi, Real* d2Ndthetadzeta, Real* d2Ndxidzeta);

			static void evalLaplacian(const Real* xi, Real* lapN);

		}; // class LagrangeQuad

	} // namespace fem

} // namespace pdesolver

#endif
