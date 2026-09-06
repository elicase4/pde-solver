#ifndef PDESOLVER_FEM_BOUNDARY_NATURALBOUNDARYREGISTRY_HPP
#define PDESOLVER_FEM_BOUNDARY_NATURALBOUNDARYREGISTRY_HPP

#include <memory>
#include <unordered_map>
#include <vector>

#include "BoundaryCondition.hpp"
#include "NaturalBoundaryOperator.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			template<typename EvalQP>
			class NaturalBoundaryRegistry {

			public:

				template<BoundaryFunction Function, typename FormRegistry, typename Model>
				void registerBC(std::shared_ptr<BoundaryCondition<Function>> bc, FormRegistry& forms, Model& model) {
					entries_[bc->tag].push_back(std::make_unique<NaturalBoundaryOperator<EvalQP, Function, FormRegistry, Model>>(std::move(bc), forms, model));
				}

				void apply(Int tag, EvalQP& qp, Real* Fe) const {

					const auto* entries = getEntries(tag);

					if (!entries)
						return;

					for (const auto& entry : *entries)
						entry->apply(qp, Fe);
				}

				// per-face, before the quadrature-point loop -- mirrors apply()'s entry iteration.
				// Default no-op per entry (see NaturalBoundaryOperatorBase::gatherFaceElementData).
				void gatherFaceElementData(Int tag, const Index* faceNodeGlobalIDs, Index nodesPerFace) const {

					const auto* entries = getEntries(tag);

					if (!entries)
						return;

					for (const auto& entry : *entries)
						entry->gatherFaceElementData(faceNodeGlobalIDs, nodesPerFace);
				}

				bool hasAny(Int tag) const {
					return entries_.contains(tag);
				}

				const std::vector<std::unique_ptr<NaturalBoundaryOperatorBase<EvalQP>>>* getEntries(Int tag) const {

					auto it = entries_.find(tag);

					if (it == entries_.end())
						return nullptr;

					return &it->second;
				}
				
				bool isNatural(Int tag, Index component) const {
					
					const auto* entries = getEntries(tag);
					
					if (!entries) {
						return false;
					}

					for (const auto& entry : *entries) {
						if (entry->componentType(component) == BCCategory::Natural){
							return true;
						}
					}

					return false;
				
				}


			private:

				std::unordered_map<Int, std::vector<std::unique_ptr<NaturalBoundaryOperatorBase<EvalQP>>>> entries_;

			}; // class NaturalBoundaryRegistry

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
