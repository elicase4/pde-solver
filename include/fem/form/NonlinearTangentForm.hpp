#ifndef PDESOLVER_NONLINEARTANGENTFORM_HPP
#define PDESOLVER_NONLINEARTANGENTFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<typename Form, Int Dim, Int NodesPerElement>
			concept NonlinearTangentForm = requires (const fem::eval::ElementEval<Dim, NodesPerElement>& eleEval, Real* Ke, const Real* Ue, Real* Oe) {
				{ Form::computeElementTangentMatrix<NodesPerElement>(eleEval, Ue, Ke) } -> std::same_as<void>;
				{ Form::computeElementTangentOperator<NodesPerElement>(eleEval, Ue, Oe) } -> std::same_as<void>;
			}; // concept NonlinearTangentForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
