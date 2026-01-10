#ifndef PDESOLVER_LINEARFORM_HPP
#define PDESOLVER_LINEARFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, Int Dim, Int NodesPerElement>
			concept LinearForm =
			requires (
				const fem::eval::ElementEval<Dim, NodesPerElement>& eleEval,
				Real* Fe
			) {
				{ Form::template computeElementVector<NodesPerElement>(eleEval, Fe) };
			}; // concept LinearForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
