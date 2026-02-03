#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"
#include "fem/geometry/JacobianTransform.hpp"

namespace pdesolver::fem::form {

	template<Int SpatialDim>
	struct PoissonBilinearForm<SpatialDim> {

		template<typename EvalContext>
		PDE_HOST PDE_DEVICE static void computeElementMatrix(const EvalContext& ctx, Real* Ke){
			
			Real integrand;
			
			// element matrix assembly contribution
			for (Index a = 0; a < ctx.NumNodes; ++a){
				for (Index b = 0; b < ctx.NumNodes; ++b){
					
					integrand = 0.0;
					for (Index sD = 0; sD < SpatialDim; ++sD){
						integrand += ctx.dNdx[a*SpatialDim + sD] * ctx.dNdx[b*SpatialDim + sD];
					}

					Ke[a * ctx.NumNodes + b] += integrand * ctx.measure * ctx.w;
				}
			}
		}
		
		template<typename EvalContext>
		PDE_HOST PDE_DEVICE static void computeElementOperator(const EvalContext& ctx, const Real* Ue, Real* Oe){
			
			Real integrand;

			// element matrix assembly contribution
			for (Index a = 0; a < ctx.NumNodes; ++a){
				for (Index b = 0; b < ctx.NumNodes; ++b){
					
					integrand = 0.0;
					for (Index sD = 0; sD < SpatialDim; ++sD){
						integrand += ctx.dNdx[a*SpatialDim + sD] * ctx.dNdx[b*SpatialDim + sD];
					}
					
					Oe[a] += integrand * Ue[b] * ctx.measure * ctx.w;
				}
			}
		}

	}; // struct PoissonBilinearForm<SpatialDim>

} // namespace pdesolver::fem::form

#endif
