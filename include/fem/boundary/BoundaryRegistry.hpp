#ifndef PDESOLVER_BOUNDARYREGISTRY_HPP
#define PDESOLVER_BOUNDARYREGISTRY_HPP

#include <vector>
#include <memory>
#include <variant>
#include <unordered_map>

#include "core/Types.hpp"
#include "fem/boundary/BoundaryCondition.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			template<typename... EvalQPs>
			class BoundaryRegistry {
			public:
				
				template<typename EvalQP>
				struct BoundaryEvalQP {
					const EvalQP& qp;
					Real* Ue;
					Real* Fe;
				};

				class BCEntryBase {
				public:

					virtual ~BCEntryBase() = default;

					virtual Int tag() const = 0;
					virtual Index numComponents() const = 0;
					virtual BCCategory componentType(Index c) const = 0;
					virtual void eval(Real time, const Real* x, Real* out) const = 0;
					virtual void apply(const std::variant<BoundaryEvalQP<EvalQPs>...>& ctx) const = 0;
					
				}; // class BCEntryBase
				
				template<typename Function, typename FormRegistry, typename Model>
				struct BCEntry final : BCEntryBase {
				public:

					std::shared_ptr<BoundaryCondition<Function, FormRegistry, Model>> bc;

					explicit BCEntry(std::shared_ptr<BoundaryCondition<Function, FormRegistry, Model>> bcIn) : bc(std::move(bcIn)) {}

					Int tag() const override { return bc->tag; }

					Index numComponents() const override { return bc->NumComponents; }

					BCCategory componentType(Index c) const override { return bc->componentType[c]; }

					void eval(Real time, const Real* x, Real* out) const override { bc->f.eval(time, x, out); }
					
					void apply(const std::variant<BoundaryEvalQP<EvalQPs>...>& ctx) const override {
						std::visit([this](const auto& typedCtx) { applyImpl(typedCtx); }, ctx);
					}

				private:

					template<typename EvalQP>
					void applyImpl(const BoundaryEvalQP<EvalQP>& ctx) const {
						if constexpr (requires { bc->model.eval(ctx.qp); bc->model.evalGradient(ctx.qp); bc->forms.computeElementLevelVector(ctx.qp, ctx.Ue, ctx.Fe); }) {
							bc->model.eval(ctx.qp);
							bc->model.evalGradient(ctx.qp);
							bc->forms.computeElementLevelVector(ctx.qp, ctx.Ue, ctx.Fe);
						}
					}

				}; // struct BCEntry

				template<typename Function, typename FormRegistry, typename Model>
				void registerBC(std::shared_ptr<BoundaryCondition<Function, FormRegistry, Model>> bc){
					entries_[bc->tag].push_back(std::make_unique<BCEntry<Function, FormRegistry, Model>>(std::move(bc)));
				}

				void applyBoundaryContributions(Int tag, const std::variant<BoundaryEvalQP<EvalQPs>...>& ctx) const {

					const auto* bc = getEntries(tag);

					if (!bc){
						return;
					}

					for (const auto& bcEntry : *bc){
						bcEntry->apply(ctx);
					}

				}

				const std::vector<std::unique_ptr<BCEntryBase>>* getEntries(Int tag) const {
				
					auto it = entries_.find(tag);

					if (it == entries_.end()) {
						return nullptr;
					}

					return &it->second;
				
				}

				bool isEssential(Int tag, Index component) const {
					
					const auto* bc = getEntries(tag);

					if (!bc) {
						return false;
					}

					for (const auto& bcEntry : *bc) {
						if (bcEntry->componentType(component) == BCCategory::Essential) {
							return true;
						}
					} 

					return false;
				}

				bool isNatural(Int tag, Index component) const {
					
					const auto* bc = getEntries(tag);

					if (!bc) {
						return false;
					}

					for (const auto& bcEntry : *bc) {
						if (bcEntry->componentType(component) == BCCategory::Natural) {
							return true;
						}
					}

					return false;
				}

				bool hasAny(Int tag) const {
					return entries_.contains(tag);
				}

			private:
				
				std::unordered_map<Int, std::vector<std::unique_ptr<BCEntryBase>>> entries_;

			}; // class BoundaryRegistry

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
