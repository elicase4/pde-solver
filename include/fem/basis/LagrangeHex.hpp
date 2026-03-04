#ifndef PDESOLVER_LAGRANGEHEX_HPP
#define PDESOLVER_LAGRANGEHEX_HPP

#include "fem/basis/Lagrange1D.hpp"

namespace pdesolver {
	namespace fem {
		namespace basis {
			
			template <Index Px, Index Py, Index Pz>
			class LagrangeHex {
			public:
				using BasisX = Lagrange1D<Px>;
				using BasisY = Lagrange1D<Py>;
				using BasisZ = Lagrange1D<Pz>;

				PDE_HOST PDE_DEVICE static void eval(const Real* xi, Real* N);
				PDE_HOST PDE_DEVICE static void evalGradient(const Real* xi, Real* dNdxi);
				PDE_HOST PDE_DEVICE static void evalHessian(const Real* xi, Real* d2Nd2xi);
				PDE_HOST PDE_DEVICE static void evalLaplacian(const Real* xi, Real* lapN);
				
				PDE_HOST PDE_DEVICE static Real getFaceTopology(const Int rngID, Index* tangentID);
				PDE_HOST PDE_DEVICE static Index nodesPerFace(const Int rngID);
				PDE_HOST PDE_DEVICE static void getFaceNodes(const Int rngID, Index* nodeIDs);
				
				// constants
				static constexpr Index NodesPerElement = (Px + 1)*(Py + 1)*(Pz + 1);
				static constexpr Index ParamtericDim = 3;

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
