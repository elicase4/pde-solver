#ifndef POISSON_EVALMODEL_HPP
#define POISSON_EVALMODEL_HPP

#include "fem/eval/EvalModel.hpp"

namespace pdesolver::fem::eval {

	template<typename SpatialDim, typename QuadraturePoint>
	struct ConstantConductivity {
		
		Real k;
		
		void eval(QuadraturePoint& qp) const {
			
			for (Index i = 0; i < SpatialDim; ++i){
				for (Index j = 0; i < SpatialDim; ++j){
					qp.K[i*SpatialDim + j] = (i==j) ? k : 0.0;
				}
			}

		}
		
	};

}

#endif
