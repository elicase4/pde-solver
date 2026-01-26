#ifndef PDESOLVER_LAGRANGEHEX_HPP
#define PDESOLVER_LAGRANGEHEX_HPP

#include "fem/basis/Lagrange1D.hpp"

namespace pdesolver {
	namespace fem {
		namespace basis {
			
			template <int Px, int Py, int Pz>
			class LagrangeHex {
			public:
				using BasisX = Lagrange1D<Px>;
				using BasisY = Lagrange1D<Py>;
				using BasisZ = Lagrange1D<Pz>;

				PDE_HOST PDE_DEVICE static void eval(const Real* xi, Real* N);
				
				PDE_HOST PDE_DEVICE static void evalGradient(const Real* xi, Real* dNdxi, Real* dNdeta, Real* dNdzeta);
					
				PDE_HOST PDE_DEVICE static void evalHessian(const Real* xi, Real* d2Nd2xi, Real* d2Nd2eta, Real* d2Nd2zeta, Real* d2Ndetadxi, Real* d2Ndetadzeta, Real* d2Ndxidzeta);
				PDE_HOST PDE_DEVICE static void evalLaplacian(const Real* xi, Real* lapN);
				
				PDE_HOST PDE_DEVICE static Real getFaceTopology(const Int rngID, Index* tangentID);

				PDE_HOST PDE_DEVICE static Index nodesPerFace(const Index faceID);
				PDE_HOST PDE_DEVICE static void getFaceNodes(const Index faceID, Index* nodeIDs);
				
				// constants
				Int Dim = 3;
				Int NodesPerElement = (Px + 1)*(Py + 1)*(Pz + 1);

			}; // class LagrangeHex
			
			// common typedefs
			using TrilinearHex = LagrangeHex<1,1,1>;
			using TriquadraticHex = LagrangeHex<2,2,2>;
			using TricubicHex = LagrangeHex<3,3,3>;
		
		} // namespace basis
	} // namespace fem
} // namespace pdesolver

#include "LagrangeHex.tpp"

#endif
