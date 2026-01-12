#ifndef PDESOLVER_BOUNDARYREGISTRY_HPP
#define PDESOLVER_BOUNDARYREGISTRY_HPP

#include <any>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "core/Types.hpp"
#include "fem/boundary/BoundaryCondition.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			class BoundaryRegistry {
			public:
				
				// registration
				template<typename BC>
				void registerBC(Int tag, BC bc) {
					bcStorage_[tag] = std::make_any<BC>(bc);
					bcCategory_[tag] = BC::category;
				}

				// Query
				bool hasBC(Int tag) const {
					return bcCategory_.find(tag) != bcCategory_.end();
				}

				BCCategory getBCCategory(Int tag) const {
					auto it = bcCategory_.find(tag);
					if (it == bcCategory_.end()) {
						throw std::runtime_error("No BC registered for tag " + std::to_string(tag));
					}
					return it->second;
				}
				
				bool isEssential(Int tag) const {
					return (hasBC(tag) &&getBCCategory(tag) == BCCategory::Essential);
				}
				
				bool isNatural(Int tag) const {
					return (hasBC(tag) &&getBCCategory(tag) == BCCategory::Natural);
				}

				template<typename BC>
				BC getBC(Int tag) const {
					auto it = bcStorage_.find(tag);
					if (it == bcStorage_.end()) {
						throw std::runtime_error("No BC registered for tag " + std::to_string(tag));
					}
					try {
						return std::any_cast<BC>(it->second);
					} catch (const std::bad_any_cast&) {
						throw std::runtime_error("BC type mismatch for tag " + std::to_string(tag));
					}
				}
				
				// get tags
				std::vector<Int> getTagsByCategory(BCCategory category) const {
					std::vector<Int> tags;
					for (const auto& [tag, cat] : bcCategory_) {
						if (cat == category) {
							tags.push_back(tag);
						}
					}
					return tags;
				}

				std::vector<Int> getEssentialTags() const {
					return getTagsByCategory(BCCategory::Essential);
				}

				std::vector<Int> getNaturalTags() const {
					return getTagsByCategory(BCCategory::Natural);
				}

				std::vector<Int> getAllTags() const{
					std::vector<Int> tags;
					for (const auto& [tag, _] : bcCategory_){
						tags.push_back(tag);
					}
					return tags;
				}

			private:
				std::unordered_map<Int, std::any> bcStorage_;
				std::unordered_map<Int, BCCategory> bcCategory_;
			
			}; // class BoundaryRegistry

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
