#ifndef PDESOLVER_BILINEARFORM_HPP
#define PDESOLVER_BILINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {

			template<typename Form, typename EvalContext>
			concept BilinearForm = requires (const EvalContext& ctx, Real* Ke, const Real* Ue, Real* Oe) {
				{ Form::computeElementMatrix(ctx, Ke) } -> std::same_as<void>; 
				{ Form::computeElementOperator(ctx, Ue, Oe) } -> std::same_as<void>; 
			}; // concept BilinearForm
				
		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
