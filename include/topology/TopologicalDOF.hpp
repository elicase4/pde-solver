#ifndef PDESOLVER_TOPOLOGICALDOF_HPP
#define PDESOLVER_TOPOLOGICALDOF_HPP

#include <vector>
#include <unordered_map>

#include "core/Types.hpp"
#include "fem/boundary/BoundaryRegistry.hpp"
#include "mesh/Mesh.hpp"

namespace pdesolver {
	namespace topology {

		class TopologicalDOF {
		public:
			TopologicalDOF(const mesh::Mesh& mesh, Index dofsPerNode);
			
			// size
			Index numGlobalDOFs() const { return numGlobalDOFs_; }
			Index numFreeDOFs() const { return numFreeDOFs_; }
			Index dofsPerNode() const { return dofsPerNode_; }
			
			// node mappings
			Index getNodeDOF(Index nodeId, Index component) const { return (nodeId * dofsPerNode_ + component); }
			
			// element mappings
			void getElementDOFs(Index elemId, Index* dofs) const;
			
			// constraints
			template<typename Element>
			void buildConstraints(const fem::boundary::BoundaryRegistry& bcRegistry);

			// constraint query
			bool isConstrained(Index topoDOF) const { return topoToAlg_[topoDOF] == -1; }
			Int getConstraintTag(Index topoDOF) const;

			// mappings
			Index toAlgebraic(Index topoDOF) const { return topoToAlg_[topoDOF]; }
			Index* getTopoToAlg() const { return TopoToAlg_.data(); }
			Index toTopological(Index algDOF) const { return algToTopo_[algDOF]; }
			Index* getAlgToTopo() const { return AlgToTopo_.data(); }

		private:
			const mesh::Mesh& mesh_;
			Index dofsPerNode_;
			Index numGlobalDOFs_;
			Index numFreeDOFs_;

			// mappings
			std::vector<Index> topoToAlg_; // Size: numGlobalDOFs
			std::vector<Index> algToTopo_; // Size: numFreeDOFs

			// constraints
			std::unordered_map<Index, Int> constraintTags_;

		}; // class TopologicalDOF

	} // namespace topology
} // namespace pdesolver

#include "TopologicalDOF.tpp"

#endif
