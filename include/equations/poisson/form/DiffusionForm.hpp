#ifndef POISSON_DIFFUSIONFORM_HPP
#define POISSON_DIFFUSIONFORM_HPP

#include "fem/form/BilinearForm.hpp"

namespace pdesolver::fem::form {

	template<typename QuadraturePointVolume>
	struct PoissonDiffusionForm {

		PDE_HOST PDE_DEVICE static void computeElementLevelMatrix(const QuadraturePointVolume& qp, const Real*, Real* Ke){
			
			Real integrand, matvecprod;
			
			// qpent matrix assembly contribution
			for (Index a = 0; a < qp.NodesPerElement; ++a){
				for (Index b = 0; b < qp.NodesPerElement; ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < qp.SpatialDim; ++sDi){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < qp.SpatialDim; ++sDj){
							matvecprod += qp.K[sDi*qp.SpatialDim + sDj] * qp.dNdx[b*qp.SpatialDim + sDj];
						}
						integrand += qp.dNdx[a*qp.SpatialDim + sDi] * matvecprod;
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
					for (Index sDi = 0; sDi < qp.SpatialDim; ++sDi){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < qp.SpatialDim; ++sDj){
							matvecprod += qp.K[sDi*qp.SpatialDim + sDj] * qp.dNdx[b*qp.SpatialDim + sDj];
						}
						integrand += qp.dNdx[a*qp.SpatialDim + sDi] * matvecprod;
					}
					
					Oe[a] += integrand * Ue[b] * qp.measure * qp.w;
				}
			}
		}

	}; // struct PoissonDiffusionForm

} // namespace pdesolver::fem::form

#endif
