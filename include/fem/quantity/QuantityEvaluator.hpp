#ifndef PDESOLVER_FEM_QUANTITY_QUANTITYEVALUATOR_HPP
#define PDESOLVER_FEM_QUANTITY_QUANTITYEVALUATOR_HPP

#include "config/Platform.hpp"

#include "core/Types.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/eval/EvalElement.hpp"
#include "fem/quantity/BoundaryQuantityRegistry.hpp"
#include "fem/quantity/Reduction.hpp"
#include "fem/quantity/QuantityForms.hpp"

#include "mesh/Mesh.hpp"

#include "linalg/types/Vector.hpp"

#include "topology/TopologicalDOF.hpp"

namespace pdesolver {
	namespace fem {
		namespace quantity {

			template<typename Backend>
			class QuantityEvaluator {
			public:

				// reduces forms over the whole mesh domain
				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename QuantityFormsT, typename Quadrature>
				static void evaluateDomain(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const Model& model, const QuantityFormsT& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, Backend>& U, Real* out);

				// reduces forms over combo boundary faces
				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename QuantityFormsT, typename Quadrature>
				static void evaluateBoundaryRegistry(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry& bcRegistry, const Real time, const Model& model, const QuantityFormsT& forms, const EvalEle& evalEle, const Quadrature& quadrature, const linalg::types::Vector<Real, Backend>& U, BoundaryQuantityRegistry<QuantityFormsT>& registry);

			}; // class QuantityEvaluator

		} // namespace quantity
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/QuantityEvaluator.tpp"

#endif
