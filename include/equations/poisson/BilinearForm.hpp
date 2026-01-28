#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"
#include "fem/geometry/JacobianTransform.hpp"

namespace pdesolver::fem::form {

	template<Int Dim>
	struct PoissonBilinearForm; // struct PoissonBilinearForm

	template<>
	struct PoissonBilinearForm<2> {

		template<typename EvalContext>
		PDE_HOST PDE_DEVICE static void computeElementMatrix(const EvalContext& ctx, Real* Ke){
			
			// element matrix assembly contribution
			for (Index a = 0; a < ctx.NumNodes; ++a){
				for (Index b = 0; b < ctx.NumNodes; ++b){
					Ke[a * ctx.NumNodes + b] += (ctx.dNdx[a]*ctx.dNdx[b] + ctx.dNdy[a]*ctx.dNdy[b]) * ctx.detJ * ctx.w;
				}
			}
		}
		
		template<typename EvalContext>
		PDE_HOST PDE_DEVICE static void computeElementOperator(const EvalContext& ctx, const Real* Ue, Real* Oe){
			
			// element matrix assembly contribution
			for (Index a = 0; a < ctx.NumNodes; ++a){
				for (Index b = 0; b < ctx.NumNodes; ++b){
					Oe[a] += (ctx.dNdx[a]*ctx.dNdx[b] + ctx.dNdy[a]*ctx.dNdy[b]) * Ue[b] * ctx.detJ * ctx.w;
				}
			}
		}

	}; // struct PoissonBilinearForm<2>

} // namespace pdesolver::fem::form

#endif
