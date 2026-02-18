#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "equations/poisson/PoissonSourceTerm.hpp"

#include "fem/core/Types.hpp"
#include "fem/config/Platform.hpp"
#include "fem/eval/EvalQuadraturePoint.hpp"
#include "fem/geometry/JacobianTransform.hpp"

namespace pdesolver::fem::form {

	template<Int SpatialDim, typename SourceFunction>
	requires fem::eval::EvalFunction<SourceFunction>
	struct PoissonLinearForm<SpatialDim> {
		
		SourceFunction source;

		template<typename QuadraturePoint>
		PDE_HOST PDE_DEVICE static void computeElementVector(const auto& qp, Real* Fe){
			
			Real val[1];
			source.eval(qp.time, qp.x, val);

			// element vector assembly contribution
			for (Index a = 0; a < qp.NumNodes; ++a){
				Fe[a] += (val[0] * qp.N[a]) * qp.measure * qp.w;
			}
		}

	};// struct PoissonLinearForm<SpatialDim>

} // namespace pdesolver::fem::form
