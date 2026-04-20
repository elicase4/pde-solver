#ifndef PDESOLVER_NONLINEARTANGENTFORM_HPP
#define PDESOLVER_NONLINEARTANGENTFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/EvalQuadraturePointVolume.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, typename QuadraturePoint, typename refPtr>
			concept NonlinearTangentForm = requires (const Form f, const QuadraturePoint& qp, refPtr Ue, Real* Ke, Real* Oe) {
				{ f.computeElementLevelMatrix(qp, Ue, Ke) } -> std::same_as<void>;
				{ f.computeElementLevelVector(qp, Ue, Oe) } -> std::same_as<void>;
			}; // concept NonlinearTangentForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
