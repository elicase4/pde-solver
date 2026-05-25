#ifndef PDESOLVER_MESH_EXCHANGE_GMSH_MESHCONVERTER_HPP
#define PDESOLVER_MESH_EXCHANGE_GMSH_MESHCONVERTER_HPP

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/exchange/gmsh/IntermediateMesh.hpp"
#include "mesh/exchange/gmsh/Utils.hpp"

namespace pdesolver {
	namespace mesh {
		namespace exchange {
			namespace gmsh {

				class MeshConverter {
				public:
					
					static void toSolverMesh(mesh::Mesh& mesh, const IntermediateMesh& input, const std::unordered_map<Int, Int>& physicalGroupMap = {});

				private:

					static ElementType buildConnectivity(mesh::Mesh& mesh, const IntermediateMesh& input);

					static void buildBoundaryTags(mesh::Mesh& mesh, const IntermediateMesh& input, ElementType cellType, const std::unordered_map<Int, Int>& physicalGroupMap);

					static std::vector<Index> reorderConnectivity(const Index* conn, ElementType type);

					static std::vector<Index> localFaceNodesSolver(const Index* elemNodes, ElementType type, Index face);

				}; // class MeshConverter

			} // namespace gmsh
		} // namespace exchange
	} // namespace mesh
} // namespace pdesolver

#endif
