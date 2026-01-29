#ifndef PDESOLVER_NONLINEARTANGENTFORM_HPP
#define PDESOLVER_NONLINEARTANGENTFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, typename EvalContext>
			concept NonlinearTangentForm = requires (const EvalContext& ctx, Real* Ke, const Real* Ue, Real* Oe) {
				{ Form::computeElementMatrix(ctx, Ue, Ke) } -> std::same_as<void>;
				{ Form::computeElementOperator(ctx, Ue, Oe) } -> std::same_as<void>;
			}; // concept NonlinearTangentForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
