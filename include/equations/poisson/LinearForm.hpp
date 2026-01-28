#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "equations/poisson/PoissonSourceTerm.hpp"

#include "fem/core/Types.hpp"
#include "fem/config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"
#include "fem/geometry/JacobianTransform.hpp"

namespace pdesolver::fem::form {

	template<Int Dim>
	struct PoissonLinearForm; // struct PoissonBilinearForm
	
	template<>
	struct PoissonLinearForm<2> {

		template<typename EvalContext>
		PDE_HOST PDE_DEVICE static void computeElementVector(const EvalContext& cxt, Real* Fe){
			
			// element vector assembly contribution
			for (Index a = 0; a < ctx.NumNodes; ++a){
				Fe[a] += (fem::eval::PoissonSourceTerm<2>::value(ctx.t, ctx.x) * ctx.N[a]) * ctx.detJ * ctx.w;
			}
		}

	}; // struct PoissonLinearForm<2>

} // namespace pdesolver::fem::form
