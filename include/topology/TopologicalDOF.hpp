#ifndef PDESOLVER_TOPOLOGICALDOF_HPP
#define PDESOLVER_TOPOLOGICALDOF_HPP

#include "mesh/Mesh.hpp"
#include "core/Types.hpp"

namespace pdesolver {
	namespace fem {
		namespace dof {

			class TopologicalDOF {
			public:
				TopologicalDOF(const mesh::Mesh& mesh, Index dofsPerNode) : mesh_(mesh), dofsPerNode_(dofsPerNode), numGlobalDOFs_ {};
				
				// size
				Index numGlobalDOFs() const { return numGlobalDOFs_; }
				Index dofsPerNode() const { return dofsPerNode_; }
				
				// node mappings
				Index getNodeDOF(Index nodeId, Index component) const { return (nodeId * dofsPerNode_ + component); }
				Index* getNodeDOFs(Index nodeId) { return (nodeId * dofsPerNode_); }
				
				// element mappings
				void getElementDOFs(Index elemId, Index* dofs) const {
					const Index* nodes = mesh_.getElementNodes(elemId);
					const Index npe = mesh_.data.nodesPerElement;

					Index k = 0;
					for (Index a = 0; a < npe; ++a) {
						for (Index c = 0; c < dofsPerNode_; ++c) {
							dofs[k++] = nodes[a] * dofsPerNode_ + c;
						}
					}
				}

			private:
				const mesh::Mesh& mesh_;
				Index dofsPerNode_;
				Index numGlobalDOFs_;

			}; // class TopologicalDOF

		} // namespace dof
	} // namespace fem
} // namespace pdesolver

#endif
