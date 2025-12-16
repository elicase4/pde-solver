#ifndef PDESOLVER_LAGRANGEQUAD_HPP
#define PDESOLVER_LAGRANGEQUAD_HPP

#include "fem/Lagrange1D.hpp"

namespace pdesolver {
	
	namespace fem {
		
		template <int Px, int Py>
		class LagrangeQuad {
		public:
			static constexpr int nodesPerElement = (Px + 1) * (Py + 1);

			static void eval(const Real* xi, Real* N);

			static void evalGradient(const Real* xi, Real* dNdxi, Real* dNdtheta);

			static void evalHessian(const Real* xi, Real* d2Nd2xi, Real* d2Nd2theta, Real* d2Ndthetadxi);

			static void evalLaplacian(const Real* xi, Real* lapN);

		}; // class LagrangeQuad

	} // namespace fem

} // namespace pdesolver

#endif
