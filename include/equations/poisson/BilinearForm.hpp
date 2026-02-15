#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"
#include "fem/geometry/JacobianTransform.hpp"

namespace pdesolver::fem::form {

	template<Int SpatialDim>
	struct PoissonBilinearForm<SpatialDim> {

		template<typename EvalElement>
		PDE_HOST PDE_DEVICE static void computeElementMatrix(const EvalContext& ctx, Real* Ke){
			
			Real integrand, matvecprod;
			
			// element matrix assembly contribution
			for (Index a = 0; a < ctx.NumNodes; ++a){
				for (Index b = 0; b < ctx.NumNodes; ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < SpatialDim; ++sD){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < SpatialDim; ++sD){
							matvecprod += ctx.K[dSi*SpatialDim + sDj] * ctx.dNdx[b*SpatialDim + sDj]
						}
						integrand += ctx.dNdx[a*SpatialDim + sDi] * innerprod;
					}

					Ke[a * ctx.NumNodes + b] += integrand * ctx.measure * ctx.w;
				}
			}
		}
		
		template<typename EvalContext>
		PDE_HOST PDE_DEVICE static void computeElementOperator(const EvalContext& ctx, const Real* Ue, Real* Oe){
			
			Real integrand, matvecprod;

			// element matrix assembly contribution
			for (Index a = 0; a < ctx.NumNodes; ++a){
				for (Index b = 0; b < ctx.NumNodes; ++b){
					
					integrand = 0.0;
					for (Index sDi = 0; sDi < SpatialDim; ++sD){
						matvecprod = 0.0;
						for (Index sDj = 0; sDj < SpatialDim; ++sD){
							matvecprod += ctx.K[dSi*SpatialDim + sDj] * ctx.dNdx[b*SpatialDim + sDj]
						}
						integrand += ctx.dNdx[a*SpatialDim + sDi] * innerprod;
					}
					
					Oe[a] += integrand * Ue[b] * ctx.measure * ctx.w;
				}
			}
		}

	}; // struct PoissonBilinearForm<SpatialDim>

} // namespace pdesolver::fem::form

#endif
