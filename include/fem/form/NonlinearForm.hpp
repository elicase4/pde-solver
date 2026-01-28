#ifndef PDESOLVER_NONLINEARFORM_HPP
#define PDESOLVER_NONLINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, typename EvalContext>
			concept NonlinearForm = requires (const EvalContext& ctx, const Real* Ue, Real* Re) {
				{ Form::computeElementResidual(ctx, Ue, Re) } -> std::same_as<void>;
			}; // concept NonlinearForm
		
		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
