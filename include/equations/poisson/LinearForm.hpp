#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "equations/poisson/PoissonSourceTerm.hpp"

#include "fem/core/Types.hpp"
#include "fem/config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"
#include "fem/geometry/JacobianTransform.hpp"

namespace pdesolver::fem::form {

	template<>
	struct PoissonLinearForm<SpatialDim> {

		template<typename EvalContext>
		PDE_HOST PDE_DEVICE static void computeElementVector(const EvalContext& ctx, Real* Fe){
			
			// element vector assembly contribution
			for (Index a = 0; a < ctx.NumNodes; ++a){
				Fe[a] += (ctx.rhsF * ctx.N[a]) * ctx.measure * ctx.w;
			}
		}

	};// struct PoissonLinearForm<SpatialDim>

} // namespace pdesolver::fem::form
