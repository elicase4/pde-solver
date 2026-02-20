#ifndef PDESOLVER_BILINEARFORM_HPP
#define PDESOLVER_BILINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/EvalQuadraturePoint.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {

			template<typename Form, typename QuadraturePoint>
			concept BilinearForm = requires (const QuadraturePoint& qp, Real* Ke, const Real* Ue, Real* Oe) {
				{ Form::computeElementMatrix(qp, Ke) } -> std::same_as<void>; 
				{ Form::computeElementOperator(qp, Ue, Oe) } -> std::same_as<void>; 
			}; // concept BilinearForm
				
		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
