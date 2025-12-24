#ifndef PDESOLVER_LINEARFORM_HPP
#define PDESOLVER_LINEARFORM_HPP

#include "core/Types.hpp"
#include "core/CudaMacros.hpp"
#include "fem/eval/ElementEval.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {
			
			template<Int Dim, Int NodesPerElement>
			class LinearForm {
				HOST_DEVICE static void computeElementVector(const eval::ElementEval<Dim,NodesPerElement>& eleEval, Real* fe);
			}; // class LinearForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
