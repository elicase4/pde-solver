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

				PDE_HOST PDE_DEVICE static void eval(const Real* xi, Real* N);
				PDE_HOST PDE_DEVICE static void evalGradient(const Real* xi, Real* dNdxi);
				PDE_HOST PDE_DEVICE static void evalHessian(const Real* xi, Real* d2Nd2xi);
				PDE_HOST PDE_DEVICE static void evalLaplacian(const Real* xi, Real* lapN);
				
				PDE_HOST PDE_DEVICE static Real getFaceTopology(const Int faceID, Index* tangentID);
				PDE_HOST PDE_DEVICE static Index nodesPerFace(const Index faceID);
				PDE_HOST PDE_DEVICE static void getFaceNodes(const Index faceID, Index* nodeIDs);
				
				// constants
				Int NodesPerElement = (Px + 1)*(Py + 1);
				Int ParametricDim = 2;
			
			}; // class LagrangeQuad

			// common typedefs
			using BilinearQuad = LagrangeQuad<1,1>;
			using BiquadraticQuad = LagrangeQuad<2,2>;
			using BicubicQuad = LagrangeQuad<3,3>;
			
		} // namespace basis
	} // namespace fem
} // namespace pdesolver

#include "LagrangeQuad.tpp"

#endif
