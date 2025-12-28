#ifndef POISSON_SOURCETERM_HPP
#define POISSON_SOURCETERM_HPP

#include <cmath>

#include "fem/eval/FieldEval.hpp"

using namespace pdesolver;

template<Int Dim>
class PoissonSourceTerm: public fem::eval::ScalarFieldEval<Dim> {};

template<>
class PoissonSourceTerm<2> {
	HOST_DEVICE Real value(const Real, const Real* x){
		return (2 * M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]));
	}
};

#endif
