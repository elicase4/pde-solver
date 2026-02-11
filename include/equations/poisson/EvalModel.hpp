#ifndef POISSON_EVALMODEL_HPP
#define POISSON_EVALMODEL_HPP

#include "fem/eval/EvalContext.hpp"

namespace pdesolver::fem::eval {

	template<Int SpatialDim>
	struct ConstantConductivity {
		
		Real k;

		void evaluate(const EvalContext& ctx) const {
			
			for (Index i = 0; i < SpatialDim; ++i){
				ctx.K[i*SpatialDim + i] = k;
			}

		}
		
		void derivative(const EvalContext& ctx) const {
			
			for (Index i = 0; i < SpatialDim; ++i){
				ctx.dK[i*SpatialDim + i] = 0.0;
			}

		}

	};

}

#endif
