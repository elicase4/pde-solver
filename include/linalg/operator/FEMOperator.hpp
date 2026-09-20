#ifndef PDESOLVER_LINALG_FEMOPERATOR_HPP
#define PDESOLVER_LINALG_FEMOPERATOR_HPP

#include <array>

#include "fem/assembly/Assembler.hpp"
#include "fem/assembly/ElementMap.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"

#include "fem/eval/EvalElement.hpp"
#include "fem/eval/EvalQuadraturePointVolume.hpp"
#include "fem/eval/EvalModel.hpp"

#include "fem/form/BilinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"

#include "mesh/Mesh.hpp"

#include "linalg/types/Vector.hpp"

#include "topology/TopologicalDOF.hpp"

namespace pdesolver {
	namespace linalg {
		namespace op {

			template<typename Assembler, typename TopologicalDOF, fem::eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature, fem::assembly::GatherMode Mode, typename VectorType>
			class FEMOperator {
			public:

				static constexpr Index NumAuxStates = EvalQP::NumAuxStates;

				FEMOperator(const Assembler& assembler_, const mesh::Mesh& mesh_, const TopologicalDOF& topoDOF_, const Real* time_, const Model& model_, const FormRegistry& forms_, const EvalEle& evalEle_, const Quadrature& quadrature_, const fem::boundary::EssentialBoundaryRegistry* bcRegistry_, const VectorType* fieldSource_ = nullptr, const std::array<const VectorType*, NumAuxStates>& auxStates_ = {}) : assembler(assembler_), mesh(mesh_), topoDOF(topoDOF_), time(time_), model(model_), forms(forms_), quadrature(quadrature_), evalEle(evalEle_), bcRegistry(bcRegistry_), fieldSource(fieldSource_), auxStates(auxStates_) {}

				void apply(const VectorType& x, VectorType& y) const {
					assembler.template assembleVector<TopologicalDOF::dofsPerNode, EvalEle, EvalQP, Model, FormRegistry, Quadrature, Mode>(mesh, topoDOF, *time, model, forms, evalEle, quadrature, x, fieldSource, auxStates, y, bcRegistry);
				}

				Index size() const {
					return topoDOF.numFreeDOFs();
				}

				// Approximation by dominant local dense-matvec term only
				Index flopsPerApply() const {
					return 2 * mesh.data.numElements * mesh.data.nodesPerElement * mesh.data.nodesPerElement;
				}

			private:
				const Assembler& assembler;
				const mesh::Mesh& mesh;
				const TopologicalDOF& topoDOF;
				const Real* time;
				const Model& model;
				const FormRegistry& forms;
				const Quadrature& quadrature;
				const EvalEle& evalEle;

				const fem::boundary::EssentialBoundaryRegistry* bcRegistry;
				const VectorType* fieldSource;
				std::array<const VectorType*, NumAuxStates> auxStates;

			}; // class FEMOperator

		} // namespace op
	} // namespace linalg
} // namespace pdesolver

#endif
