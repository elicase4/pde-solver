#ifndef POISSON_DEFAULTMODEL_HPP
#define POISSON_DEFAULTMODEL_HPP

#include "fem/eval/EvalModel.hpp"

namespace pdesolver::fem::eval {

	template<typename QuadraturePoint>
	struct PoissonDefaultModel {
		
		void eval(QuadraturePoint&) const {}

		void evalGradient(QuadraturePoint&) const {}
		
	}; // struct PoissonDefaultModel

}

#endif
