#ifndef POISSON_DEFAULTMODEL_HPP
#define POISSON_DEFAULTMODEL_HPP

#include "fem/eval/EvalModel.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePoint>
	struct DefaultModel {
		
		void eval(QuadraturePoint&) const {}

		void evalGradient(QuadraturePoint&) const {}
		
	}; // struct DefaultModel

} // pdesolver::equations::heateq

#endif
