#ifndef PDESOLVER_FIELDEVAL_HPP
#define PDESOLVER_FIELDEVAL_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<typename FieldEval, Int Dim>
			concept ScalarFieldEval =
			requires (
				const Real t,
				const Real* x
			) {
				{ FielEval::template value(t, x) };
			}; // concept ScalarFieldEval

			template<typename FieldEval, Int Dim>
			concept VectorFieldEval =
			requires (
				const Real t,
				const Real* x,
				Real* v
			) {
				{ FieldEval::template value(t, x, v) };
			}; // concept VectorFieldEval
			
		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
