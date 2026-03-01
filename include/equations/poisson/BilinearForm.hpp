#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "fem/form/BilinearForm.hpp"

namespace pdesolver::fem::form {

	template<Int SpatialDim>
	struct PoissonBilinearForm {

		template<typename QuadraturePoint>
		PDE_HOST PDE_DEVICE static void computeElementMatrix(const auto& qp, Real* Ke){
			
			Real integrand, matvecprod;
			
			// qpent matrix assembly contribution
			for (Index a = 0; a < qp.NumNodes; ++a){
				for (Index b = 0; b < qp.NumNodes; ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < SpatialDim; ++sDi){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < SpatialDim; ++sDj){
							matvecprod += qp.K[sDi*SpatialDim + sDj] * qp.dNdx[b*SpatialDim + sDj];
						}
						integrand += qp.dNdx[a*SpatialDim + sDi] * matvecprod;
					}

					Ke[a * qp.NumNodes + b] += integrand * qp.measure * qp.w;
				}
			}
		}
		
		template<typename QuadraturePoint>
		PDE_HOST PDE_DEVICE static void computeElementOperator(const auto& qp, const Real* Ue, Real* Oe){
			
			Real integrand, matvecprod;

			// qpent matrix assembly contribution
			for (Index a = 0; a < qp.NumNodes; ++a){
				for (Index b = 0; b < qp.NumNodes; ++b){
					
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

	}; // struct PoissonBilinearForm<SpatialDim>

} // namespace pdesolver::fem::form

#endif
