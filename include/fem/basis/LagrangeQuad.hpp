#ifndef PDESOLVER_LAGRANGEQUAD_HPP
#define PDESOLVER_LAGRANGEQUAD_HPP

#include "fem/basis/Lagrange1D.hpp"

namespace pdesolver {
	namespace fem {
		namespace basis {

			template <int Px, int Py>
			class LagrangeQuad {
			public:
				using BasisX = Lagrange1D<Px>;
				using BasisY = Lagrange1D<Py>;

				HOST_DEVICE static void eval(const Real* xi, Real* N);
				HOST_DEVICE static void evalGradient(const Real* xi, Real* dNdxi, Real* dNdeta);
				HOST_DEVICE static void evalHessian(const Real* xi, Real* d2Nd2xi, Real* d2Nd2eta, Real* d2Ndetadxi);
				HOST_DEVICE static void evalLaplacian(const Real* xi, Real* lapN);

			}; // class LagrangeQuad
			
			// Implementation: eval
			template<int Px, int Py>
			HOST_DEVICE void LagrangeQuad<Px, Py>::eval(const Real* xi, Real* N){
				Real Nx[Px + 1];
				Real Ny[Py + 1];
				
				BasisX::eval(xi[0], Nx);
				BasisY::eval(xi[1], Ny);

				// compute tensor product
				Index a = 0;
				for (Index j = 0; j <= Py; ++j){
					for (Index i = 0; i <= Px; ++i){
						N[a] = Nx[i] * Ny[j];
						a++;
					}
				}
			}

			// Implementation: evalGradient
			template<int Px, int Py>
			HOST_DEVICE void LagrangeQuad<Px, Py>::evalGradient(const Real* xi, Real* dNdxi, Real* dNdeta){
				Real Nx[Px + 1], Ny[Py + 1];
				Real dNx[Px + 1], dNy[Py + 1];
				
				BasisX::eval(xi[0], Nx);
				BasisX::evalFirstDerivative(xi[0], dNx);
				BasisY::eval(xi[1], Ny);
				BasisY::evalFirstDerivative(xi[1], dNy);

				// evaluate tensor product & chain rule
				Index a = 0;
				for (Index j = 0; j <= Py; ++j){
					for (Index i = 0; i <= Px; ++ i){
						dNdxi[a] = dNx[i] * Ny[j];
						dNdeta[a] = Nx[i] * dNy[j];
						a++;
					}
				}
			}

			// Implementation: evalHessian
			template<int Px, int Py>
			HOST_DEVICE void LagrangeQuad<Px, Py>::evalHessian(const Real* xi, Real* d2Nd2xi, Real* d2Nd2eta, Real* d2Ndetadxi){
				Real Nx[Px + 1], Ny[Py + 1];
				Real dNx[Px + 1], dNy[Py + 1];
				Real d2Nx[Px + 1], d2Ny[Py + 1];
				
				BasisX::eval(xi[0], Nx);
				BasisX::evalFirstDerivative(xi[0], dNx);
				BasisX::evalSecondDerivative(xi[0], d2Nx);
				BasisY::eval(xi[1], Ny);
				BasisX::evalFirstDerivative(xi[1], dNy);
				BasisY::evalSecondDerivative(xi[1], d2Ny);

				// evaluate tensor product & chain rule
				Index a = 0;
				for (Index j = 0; j <= Py; ++j){
					for (Index i = 0; i <= Px; ++ i){
						d2Nd2xi[a] = d2Nx[i] * Ny[j];
						d2Nd2eta[a] = Nx[i] * d2Ny[j];
						d2Ndetadxi[a] = dNx[i] * dNy[j];
						a++;
					}
				}
			}
			
			// Implementation: evalLaplacian
			template<int Px, int Py>
			HOST_DEVICE void LagrangeQuad<Px, Py>::evalLaplacian(const Real* xi, Real* lapN){
				Real Nx[Px + 1], Ny[Py + 1];
				Real d2Nx[Px + 1], d2Ny[Py + 1];
				
				BasisX::eval(xi[0], Nx);
				BasisX::evalSecondDerivative(xi[0], d2Nx);
				BasisY::eval(xi[1], Ny);
				BasisY::evalSecondDerivative(xi[1], d2Ny);

				// evaluate tensor product & chain rule
				Index a = 0;
				for (Index j = 0; j <= Py; ++j){
					for (Index i = 0; i <= Px; ++ i){
						lapN[a] = d2Nx[i] * Ny[j] + Nx[i] * d2Ny[j];
						a++;
					}
				}
			}
			
			// common typedefs
			using BilinearQuad = LagrangeQuad<1,1>;
			using BiquadraticQuad = LagrangeQuad<2,2>;
			using BicubicQuad = LagrangeQuad<3,3>;

		} // namespace basis
	} // namespace fem
} // namespace pdesolver

#endif
