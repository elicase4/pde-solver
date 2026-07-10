#ifndef PDESOLVER_BOUNDARYREGISTRY_HPP
#define PDESOLVER_BOUNDARYREGISTRY_HPP

#include <vector>
#include <memory>

#include "core/Types.hpp"
#include "fem/boundary/BoundaryCondition.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			class BoundaryRegistry {
			public:
				
				template<typename Function, typename FormRegistry, typename Model>
				struct BCEntry {

					BoundaryCondition<Function> bc;

					BCEntry(const BoundaryCondition<Function>& bcIn) : bc(bcIn) {}

					Int tag() const {
						return bc.tag;
					}

					Index numComponents() const {
						return bc.NumComponents;
					}

					BCCategory componentType(Index c) const {
						return bc.componentType[c];
					}

					void eval(Real time, const Real* x, Real* out) const {
						bc.f.eval(time, x, out);
					}
					
					template<typename EvalQP>
					void evalElementLevelVector(EvalQP qp, const Real* Ue, Real* Fe) {
						bc.model->eval(qp);
						bc.model->evalGradient(qp);
						bc.forms->computeElementLevelVector(qp, Ue, Fe);
					}

				}; // struct BCEntry

				template<typename Function, typename FormRegistry, typename Model>
				void registerBC(const BoundaryCondition<Function>& bc){
					entries_.push_back(std::make_unique<BCEntry<Function, FormRegistry, Model>>(bc));
				}

				const auto& entries() const {
					return entries_;
				}

				bool isEssential(Int tag, Index component) const {
					
					for (const auto& bc: entries_){
						if (bc->tag() != tag) continue;
						if (bc->componentType(component) == BCCategory::Essential) {
							return true;
						}
					}

					return false;
				}

				bool isNatural(Int tag, Index component) const {
					
					for (const auto& bc: entries_){
						if (bc->tag() != tag) continue;
						if (bc->componentType(component) == BCCategory::Natural) {
							return true;
						}
					}

					return false;
				}

				bool hasAny(Int tag) const {
					
					for (const auto& bc : entries_) {
						if (bc->tag() == tag) return true;
					}

					return false;
				}

			private:
				
				std::vector<std::unique_ptr<BCEntryBase>> entries_;

						
			}; // class BoundaryRegistry

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
