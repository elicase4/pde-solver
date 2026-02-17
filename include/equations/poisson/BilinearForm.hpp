#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/qpentEval.hpp"
#include "fem/geometry/JacobianTransform.hpp"

namespace pdesolver::fem::form {

	template<Int SpatialDim>
	struct PoissonBilinearForm<SpatialDim> {

		template<typename Evalqpent>
		PDE_HOST PDE_DEVICE static void computeqpentMatrix(const auto& qp, Real* Ke){
			
			Real integrand, matvecprod;
			
			// qpent matrix assembly contribution
			for (Index a = 0; a < qp.NumNodes; ++a){
				for (Index b = 0; b < qp.NumNodes; ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < SpatialDim; ++sD){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < SpatialDim; ++sD){
							matvecprod += qp.K[dSi*SpatialDim + sDj] * qp.dNdx[b*SpatialDim + sDj]
						}
						integrand += qp.dNdx[a*SpatialDim + sDi] * innerprod;
					}

					Ke[a * qp.NumNodes + b] += integrand * qp.measure * qp.w;
				}
			}
		}
		
		template<typename Evalqpent>
		PDE_HOST PDE_DEVICE static void computeqpentOperator(const Evalqpent& qp, const Real* Ue, Real* Oe){
			
			Real integrand, matvecprod;

			// qpent matrix assembly contribution
			for (Index a = 0; a < qp.NumNodes; ++a){
				for (Index b = 0; b < qp.NumNodes; ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < SpatialDim; ++sD){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < SpatialDim; ++sD){
							matvecprod += qp.K[dSi*SpatialDim + sDj] * qp.dNdx[b*SpatialDim + sDj]
						}
						integrand += qp.dNdx[a*SpatialDim + sDi] * innerprod;
					}
					
					Oe[a] += integrand * Ue[b] * qp.measure * qp.w;
				}
			}
		}

	}; // struct PoissonBilinearForm<SpatialDim>

} // namespace pdesolver::fem::form

#endif
