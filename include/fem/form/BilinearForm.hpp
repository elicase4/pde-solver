#ifndef PDESOLVER_BILINEARFORM_HPP
#define PDESOLVER_BILINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {

			template<typename Form, Int Dim, Int NodesPerElement>
			concept BilinearForm = requires (const fem::eval::ElementEval<Dim, NodesPerElement>& eleEval, Real* Ke, const Real* Ue, Real* Oe) {
				{ Form::computeElementMatrix<NodesPerElement>(eleEval, Ke) } -> std::same_as<void>; 
				{ Form::computeElementOperator<NodesPerElement>(eleEval, Ue, Oe) } -> std::same_as<void>; 
			}; // concept BilinearForm
				
		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
