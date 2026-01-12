#ifndef PDESOLVER_LINEARFORM_HPP
#define PDESOLVER_LINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, Int Dim, Int NodesPerElement>
			concept LinearForm = requires (const fem::eval::ElementEval<Dim, NodesPerElement>& eleEval, Real* Fe) {
				{ Form::computeElementVector<NodesPerElement>(eleEval, Fe) } -> std::same_as<void>;
			}; // concept LinearForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
