#ifndef PDESOLVER_LAGRANGE1D_HPP
#define PDESOLVER_LAGRANGE1D_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace basis {
		
			template <Index Order>
			class Lagrange1D {
			public:
				PDE_HOST PDE_DEVICE static void eval(Real xi, Real* N);
				PDE_HOST PDE_DEVICE static void evalFirstDerivative(Real xi, Real* dN);
				PDE_HOST PDE_DEVICE static void evalSecondDerivative(Real xi, Real* d2N);
			
				// constants
				static constexpr Index NodesPerElement = (Order + 1);
				static constexpr Index ParametricDim = 1;

			}; // class Lagrange1D
			
		} // namespace basis
	} // namespace fem
} // namespace pdesolver

#include "Lagrange1D.tpp"

#endif
