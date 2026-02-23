#ifndef PDESOLVER_BOUNDARYAPPLICATOR_HPP
#define PDESOLVER_BOUNDARYAPPLICATOR_HPP

#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/boundary/BoundaryRegistry.hpp"

#include "mesh/Mesh.hpp"
#include "topology/TopologicalDOF.hpp"
#include "linalg/types/Vector.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			template<typename Backend>
			class BoundaryApplicator {
			public:

				static void apply(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& registry, linalg::types::Vector<Real, Backend>& F){
					for (Int tag: registry.getAllTags()) {
						switch (tag) {
							case (BCCategory::Essential):
								applyEssentialBC(mesh, topoDOF, registry, F);
								break;
							case (BCCategory::Natural):
								applyNaturalBC(mesh, topoDOF, registry, F);
								break;
						}
					}
				}

			
			private:

				static void applyEssentialBC(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& registry, linalg::types::Vector<Real, Backend>& F);
				
				static void applyNaturalBC(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& registry, linalg::types::Vector<Real, Backend>& F);

			}; // class BoundaryApplicator

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

//#include "backend/cpu/BoundaryApplicator.tpp"

#endif
