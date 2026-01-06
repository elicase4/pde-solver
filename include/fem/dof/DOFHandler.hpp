#ifndef PDESOLVER_DOFHANDLER_HPP
#define PDESOLVER_DOFHANDLER_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "mesh/Mesh.hpp"

#include <unordered_map>
#include <stdexcept>

namespace pdesolver {
	namespace fem  {
		namespace dof {

			class DOFHandler {
			public:
				
				// Node -> DOF indexing
				PDE_HOST PDE_DEVICE static Index getNodeDOF(const mesh::Mesh& mesh, const Index nodeId, const Index component, const Index dofsPerNode);
				PDE_HOST PDE_DEVICE static void getNodeDOFs(const mesh::Mesh& mesh, const Index nodeId, const Index dofsPerNode, Index* output);
				
				// Element -> DOF Indexing
				PDE_HOST PDE_DEVICE static void getElementDOFs(const mesh::Mesh& mesh, const Index elemId, const Index dofsPerNode, Index* output);
				
				// Boundary information
				PDE_HOST PDE_DEVICE static bool isConstrained(const Index dof, const Index* constrainedDofs, const Index numConstrained);
				PDE_HOST PDE_DEVICE static Int getConstrainedTag(const Index dof, const Index* constrainedDofs, const Index numConstrained);

			}; // class DOFHandler

		} // namespace dof
	} // namespace fem
} // namespace pdesolver

#endif
