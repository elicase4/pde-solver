#ifndef POISSON_DIFFUSIONFORM_HPP
#define POISSON_DIFFUSIONFORM_HPP

#include "fem/form/BilinearForm.hpp"

namespace pdesolver::fem::form {

	template<typename QuadraturePointVolume, Index SpatialDim>
	struct PoissonDiffusionForm {

		PDE_HOST PDE_DEVICE static void computeElementLevelMatrix(const QuadraturePointVolume& qp, const Real*, Real* Ke){
			
			Real integrand, matvecprod;
			
			// qpent matrix assembly contribution
			for (Index a = 0; a < qp.NodesPerElement; ++a){
				for (Index b = 0; b < qp.NodesPerElement; ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < SpatialDim; ++sDi){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < SpatialDim; ++sDj){
							matvecprod += qp.K[sDi*SpatialDim + sDj] * qp.dNdx[b*SpatialDim + sDj];
						}
						integrand += qp.dNdx[a*SpatialDim + sDi] * matvecprod;
					}

					Ke[a * qp.NodesPerElement + b] += integrand * qp.measure * qp.w;
				}
			}
		}
		
		PDE_HOST PDE_DEVICE static void computeElementLevelVector(const QuadraturePointVolume& qp, const Real* Ue, Real* Oe){
			
			Real integrand, matvecprod;

			// qpent matrix assembly contribution
			for (Index a = 0; a < qp.NodesPerElement; ++a){
				for (Index b = 0; b < qp.NodesPerElement; ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < SpatialDim; ++sDi){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < SpatialDim; ++sDj){
							matvecprod += qp.K[sDi*SpatialDim + sDj] * qp.dNdx[b*SpatialDim + sDj];
						}
						integrand += qp.dNdx[a*SpatialDim + sDi] * matvecprod;
					}
					
					Oe[a] += integrand * Ue[b] * qp.measure * qp.w;
				}
			}
		}

	}; // struct PoissonStiffnessMatrix

} // namespace pdesolver::fem::form

#endif
