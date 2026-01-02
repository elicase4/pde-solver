#ifndef PDESOLVER_NONLINEARTANGENTFORM_HPP
#define PDESOLVER_NONLINEARTANGENTFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<Int Dim, Int NodesPerElement>
			class NonlinearTangentForm {
				PDE_HOST PDE_DEVICE static void computeElementTangentMatrix(const eval::ElementEval<Dim,NodesPerElement>& eleEval, const Real* Ue, Real* Ke);
				PDE_HOST PDE_DEVICE static void computeElementTangentOperator(const eval::ElementEval<Dim,NodesPerElement>& eleEval, const Real* Ue, Real* Oe);
			}; // class NonlinearTangentForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
