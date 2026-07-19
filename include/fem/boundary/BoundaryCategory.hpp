#ifndef PDESOLVER_FEM_BOUNDARY_BOUNDARYCATEGORY_HPP
#define PDESOLVER_FEM_BOUNDARY_BOUNDARYCATEGORY_HPP

namespace pdesolver {
	namespace fem {
		namespace boundary {

			enum class BCCategory {
				None,
				Essential,
				Natural
			}; // enum class BCCategory

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
