#ifndef PDESOLVER_MESH_EXCHANGE_GMSH_ELEMENTTYPE_HPP
#define PDESOLVER_MESH_EXCHANGE_GMSH_ELEMENTTYPE_HPP

namespace pdesolver {
	namespace mesh {
		namespace exchange {
			namespace gmsh {

				enum class ElementType {

					Unknown = 0,

					// 1D
					LineP1,
					LineP2,

					// 2D
					TriP1,
					TriP2,

					QuadP1,
					QuadP2,

					// 3D
					TetP1,
					TetP2,

					HexP1,
					HexP2,

					// IGA
					BSplinePatch,
					NURBSPatch

				}; // enum class ElementType

			} // namespace gmsh
		} // namespace exchange
	} // namespace mesh
} // namespace pdesolver

#endif
