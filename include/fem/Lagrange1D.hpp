#ifndef PDESOLVER_LAGRANGE1D_HPP
#define PDESOLVER_LAGRANGE1D_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace fem {
		
		template <int Order>
		class Lagrange1D {
		public:
			static constexpr int numNodes = Order + 1;
		
			static void eval(Real xi, Real* N);
			static void evalFirstDerivative(Real xi, Real* dN);
			static void evalSecondDerivative(Real xi, Real* d2N);
		
		}; // class Lagrange1D

	} // namespace fem

} // namespace pdesolver

#endif
