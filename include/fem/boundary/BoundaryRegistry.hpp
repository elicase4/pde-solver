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
				
				struct BCEntryBase {
					virtual ~BCEntryBase() = default;

					virtual Int tag() const = 0;
					virtual Index numComponents() const = 0;
					
					virtual BCCategory componentType(Index c) const = 0;

					virtual void eval(Real time, const Real* x, Real* out) const = 0;
				
				}; // struct BCEntryBase

				template<typename Function>
				struct BCEntry : BCEntryBase {

					BoundaryCondition<Function> bc;

					BCEntry(const BoundaryCondition<Function>& bcIn) : bc(bcIn) {}

					Int tag() const override {
						return bc.tag;
					}

					Index numComponents() const override {
						return bc.NumComponents;
					}

					BCCategory componentType(Index c) const override {
						return bc.componentType[c];
					}

					void eval(Real time, const Real* x, Real* out) const override {
						bc.f.eval(time, x, out);
					}

				}; // struct BCEntry

				template<typename Function>
				void registerBC(const BoundaryCondition<Function>& bc){
					entries_.push_back(std::make_unique<BCEntry<Function>>(bc));
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

				template<typename Function>
				const BoundaryCondition<Function> getBC(Int tag) const {
					
					for (const auto& bc : entries_) {
						if (bc->tag() == tag) {
							return bc;
						}
					}

					return nullptr;
				}

			private:
				
				std::vector<std::unique_ptr<BCEntryBase>> entries_;

						
			}; // class BoundaryRegistry

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
