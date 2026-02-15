#ifndef POISSON_EVALMODEL_HPP
#define POISSON_EVALMODEL_HPP

#include "fem/eval/EvalContext.hpp"

namespace pdesolver::fem::eval {

	template<typename EvalElement>
	struct ConstantConductivity {
		
		Real k;

		void evaluate(const EvalElement& ctx) const {
			
			for (Index i = 0; i < EvalElement::SpatialDim; ++i){
				ctx.K[i*EvalElement::SpatialDim + i] = k;
			}

		}
		
		void derivative(const EvalElement& ctx) const {
			
			for (Index i = 0; i < EvalElement::SpatialDim; ++i){
				ctx.dK[i*EvalElement::SpatialDim + i] = 0.0;
			}

		}

	};

}

#endif
