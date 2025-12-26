#ifndef PDESOLVER_FIELDEVAL_HPP
#define PDESOLVER_FIELDEVAL_HPP

#include "core/Types.hpp"
#include "core/CudaMacros.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<Int Dim>
			class ScalarFieldEval{
				HOST_DEVICE static Real value(const Real t, const Real* x);
			}; // class ScalarFieldEval

			template<Int Dim>
			class VectorFieldEval{
				HOST_DEVICE static void value(const Real t, const Real* x, Real* v);
			}; // class VectorFieldEval

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
