#ifndef PDESOLVER_NONLINEARTANGENTFORM_HPP
#define PDESOLVER_NONLINEARTANGENTFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, Int Dim, Int NodesPerElement>
			concept NonlinearTangentForm =
			requires (
				const fem::eval::ElementEval<Dim, NodesPerElement>& eleEval,
				Real* Ke,
				const Real* Ue,
				Real* Oe
			) {
				{ Form::template computeElementTangentMatrix<NodesPerElement>(eleEval, Ue, Ke) };
				{ Form::template computeElementTangentOperator<NodesPerElement>(eleEval, Ue, Oe) };
			}; // concept NonlinearTangentForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
