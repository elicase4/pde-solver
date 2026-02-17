#ifndef PDESOLVER_NONLINEARTANGENTFORM_HPP
#define PDESOLVER_NONLINEARTANGENTFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, typename QuadraturePoint>
			concept NonlinearTangentForm = requires (const QuadraturePoint& qp, Real* Ke, const Real* Ue, Real* Oe) {
				{ Form::computeElementMatrix(qp, Ue, Ke) } -> std::same_as<void>;
				{ Form::computeElementOperator(qp, Ue, Oe) } -> std::same_as<void>;
			}; // concept NonlinearTangentForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
