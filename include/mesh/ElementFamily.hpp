#ifndef PDESOLVER_MESH_ELEMENTFAMILY_HPP
#define PDESOLVER_MESH_ELEMENTFAMILY_HPP

namespace pdesolver {
	namespace mesh {

		enum class ElementFamily {
			Quad,
			Tri,
			Hex,
			Tet,
			Wedge
		}; // enum class ElementFamily

	} // namespace mesh
} // namespace pdesolver

#endif
