#ifndef PDESOLVER_FEM_BOUNDARY_NATURALBOUNDARYOPERATOR_HPP
#define PDESOLVER_FEM_BOUNDARY_NATURALBOUNDARYOPERATOR_HPP

#include <memory>
#include <unordered_map>
#include <vector>

#include "BoundaryCondition.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			template<typename EvalQP>
			class NaturalBoundaryOperatorBase {
			public:

				virtual ~NaturalBoundaryOperatorBase() = default;
				virtual void apply(EvalQP&, Real*) const = 0;
				virtual BCCategory componentType(Index component) const = 0;

				// per-face, before the quadrature-point loop -- default no-op, so existing
				// operators need no changes. A future data-driven flux operator overrides this
				// to gather from an EvalNodalData source into its own (mutable) storage, using
				// the face-local global node IDs the caller already has in scope.
				virtual void gatherFaceElementData(const Index*, Index) const {}

			}; // class NaturalBoundaryOperatorBase

			template<typename EvalQP, BoundaryFunction Function, typename FormRegistry, typename Model>
			class NaturalBoundaryOperator final : public NaturalBoundaryOperatorBase<EvalQP> {
			public:

				NaturalBoundaryOperator(std::shared_ptr<BoundaryCondition<Function>> bc, FormRegistry& forms, Model& model) : bc_(std::move(bc)), forms_(forms), model_(model) {}

				// the flux value comes from qp.x/qp.time inside forms_ (see FluxBoundaryForm),
				// not from bc_ -- bc_ is only needed for componentType() below
				void apply(EvalQP& qp, Real* Fe) const override {
					model_.eval(qp);
					model_.evalGradient(qp);
					forms_.computeElementLevelVector(qp, nullptr, Fe);
				}

				BCCategory componentType(Index component) const override {
					return bc_->componentType[component];
				}

			private:

				std::shared_ptr<BoundaryCondition<Function>> bc_;
				FormRegistry& forms_;
				Model& model_;

			}; // class NaturalBoundaryOperator

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
