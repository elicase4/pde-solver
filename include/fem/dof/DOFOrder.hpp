#ifndef PDESOLVER_FEM_DOFORDER_HPP
#define PDESOLVER_FEM_DOFORDER_HPP

namespace pdesolver {
	namespace fem {
		namespace dof {

			enum class DOFOrdering {
				Interleaved, // node-major
				Block // field-major
			}; // enum class DOFOrdering

		} // namespace dof
	} // namespace fem
} // namespace pdesolver

#endif
