#ifndef POISSON_EVALMODEL_HPP
#define POISSON_EVALMODEL_HPP

#include "fem/eval/EvalModel.hpp"

namespace pdesolver::fem::eval {

	template<typename QuadraturePoint, Index SpatialDim>
	struct ConstantConductivity {
		
		Real k;
		
		void eval(QuadraturePoint& qp) {
			
			for (Index i = 0; i < SpatialDim; ++i){
				for (Index j = 0; j < SpatialDim; ++j){
					qp.K[i*SpatialDim + j] = (i==j) ? k : 0.0;
				}
			}

		}

		void evalGradient(QuadraturePoint& qp) {}
		
	};

}

#endif
