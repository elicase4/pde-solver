#ifndef PDESOLVER_LAGRANGEHEX_HPP
#define PDESOLVER_LAGRANGEHEX_HPP

#include "fem/Lagrange1D.hpp"

namespace pdesolver {
	
	namespace fem {
		
		template <int Px, int Py, int Pz>
		class LagrangeQuad {
		public:
			static constexpr int dim = 3;
			static constexpr int nodesPerElement = (Px + 1) * (Py + 1) * (Pz + 1);
			
			using BasisX = Lagrange1D<Px>;
			using BasisY = Lagrange1D<Py>;
			using BasisZ = Lagrange1D<Pz>;

			HOST_DEVICE static void eval(const Real* xi, Real* N);
			HOST_DEVICE static void evalGradient(const Real* xi, Real* dNdxi, Real* dNdtheta, Real* dNdzeta);
			HOST_DEVICE static void evalHessian(const Real* xi, Real* d2Nd2xi, Real* d2Nd2theta, Real* d2Nd2zeta, 
													Real* d2Ndthetadxi, Real* d2Ndthetadzeta, Real* d2Ndxidzeta);
			HOST_DEVICE static void evalLaplacian(const Real* xi, Real* lapN);

		}; // class LagrangeQuad
		
		// Implementation: eval
		template<int Px, int Py, int Pz>
		HOST_DEVICE static void LagrangeHex<Px, Py, Pz>::eval(const Real* xi, Real* N){
			Real Nx[Px + 1];
			Real Ny[Py + 1];
			Real Nz[Pz + 1];

			BasisX::eval(xi[0], Nx);
			BasisY::eval(xi[1], Ny);
			BasisZ::eval(xi[2], Nz);

			// compute tensor product
			Index a = 0;
			for (Index k = 0; k <= Pz; ++k){
				for (Index j = 0; j <= Py; ++j){
					for (Index i = 0; i <= Px; ++i){
						N[a] = Nx[i] * Ny[j] * Nz[k];
						a++;
					}
				}
			}
		}
		
		// Implementation: evalGradient
		template<int Px, int Py, int Pz>
		HOST_DEVICE static void LagrangeHex<Px, Py, Pz>::evalGradient(const Real* xi, Real* dNdxi, Real* dNdtheta, Real* dNdzeta){
			Real Nx[Px + 1], Ny[Py + 1], Nz[Pz + 1];
			Real dNx[Px + 1], dNy[Py + 1], dNz[Pz + 1];

			BasisX::eval(xi[0], Nx);
			BasisX::evalFirstDerivative(xi[0], dNx);
			BasisY::eval(xi[1], Ny);
			BasisX::evalFirstDerivative(xi[1], dNy);
			BasisZ::eval(xi[2], Nz);
			BasisX::evalFirstDerivative(xi[2], dNz);

			// compute tensor product & chain rule
			Index a = 0;
			for (Index k = 0; k <= Pz; ++k){
				for (Index j = 0; j <= Py; ++j){
					for (Index i = 0; i <= Px; ++i){
						dNdxi[a] = dNx[i] * Ny[j] * Nz[k];
						dNdtheta[a] = Nx[i] * dNy[j] * Nz[k];
						dNdzeta[a] = Nx[i] * Ny[j] * dNz[k];
						a++;
					}
				}
			}
		}
		
		// Implementation: evalHessian
		template<int Px, int Py, int Pz>
		HOST_DEVICE static void LagrangeHex<Px, Py, Pz>::evalHessian(const Real* xi, Real* d2Nd2xi, Real* d2Nd2theta, Real* d2Nd2zeta, 
																		 Real* d2Ndthetadxi, Real* d2Ndthetadzeta, Real* d2Ndxidzeta);
			Real Nx[Px + 1], Ny[Py + 1], Nz[Pz + 1];
			Real dNx[Px + 1], dNy[Py + 1], dNz[Pz + 1];
			Real d2Nx[Px + 1], d2Ny[Py + 1], d2Nz[Pz + 1];

			BasisX::eval(xi[0], Nx);
			BasisX::evalFirstDerivative(xi[0], dNx);
			BasisX::evalSecondDerivative(xi[0], d2Nx);
			BasisY::eval(xi[1], Ny);
			BasisX::evalFirstDerivative(xi[1], dNy);
			BasisX::evalSecondDerivative(xi[1], d2Ny);
			BasisZ::eval(xi[2], Nz);
			BasisX::evalFirstDerivative(xi[2], dNz);
			BasisX::evalSecondDerivative(xi[2], d2Nz);

			// compute tensor product & chain rule
			Index a = 0;
			for (Index k = 0; k <= Pz; ++k){
				for (Index j = 0; j <= Py; ++j){
					for (Index i = 0; i <= Px; ++i){
						d2Nd2xi[a] = d2Nx[i] * Ny[j] * Nz[k];
						d2Nd2theta[a] = Nx[i] * d2Ny[j] * Nz[k];
						d2Nd2zeta[a] = Nx[i] * Ny[j] * d2Nz[k];
						d2Ndthetadxi[a] = dNx[i] * dNy[j] * Nz[k];
						d2Ndxidzeta[a] = dNx[i] * Ny[j] * dNz[k];
						d2Ndthetadzeta[a] = Nx[i] * dNy[i] * dNz[k];
						a++;
					}
				}
			}
		}
		
		// Implementation: evalLaplacian
		template<int Px, int Py, int Pz>
		HOST_DEVICE static void LagrangeHex<Px, Py, Pz>::evalLaplacian(const Real* xi, Real* lapN){
			Real Nx[Px + 1], Ny[Py + 1], Nz[Pz + 1];
			Real d2Nx[Px + 1], d2Ny[Py + 1], d2Nz[Pz + 1];

			BasisX::eval(xi[0], Nx);
			BasisX::evalSecondDerivative(xi[0], d2Nx);
			BasisY::eval(xi[1], Ny);
			BasisX::evalSecondDerivative(xi[1], d2Ny);
			BasisZ::eval(xi[2], Nz);
			BasisX::evalSecondDerivative(xi[2], d2Nz);

			// compute tensor product
			Index a = 0;
			for (Index k = 0; k <= Pz; ++k){
				for (Index j = 0; j <= Py; ++j){
					for (Index i = 0; i <= Px; ++i){
						lapN[a] = d2Nx[i] * Ny[j] * Nz[k] + Nx[i] * d2Ny[j] * Nz[k] + Nx[i] * Ny[j] * d2Nz[k];
						a++;
					}
				}
			}
		}
		
		// common typedefs
		using TrilinearQuad = LangrangeHex<1,1,1>;
		using TriquadraticQuad = LangrangeHex<2,2,2>;
		using TricubicQuad = LangrangeHex<3,3,3>;

	} // namespace fem

} // namespace pdesolver

#endif
