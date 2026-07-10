#ifndef PDESOLVER_LINALG_FEMOPERATOR_HPP
#define PDESOLVER_LINALG_FEMOPERATOR_HPP

#include "fem/assembly/Assembler.hpp"

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

			template<typename Assembler, typename TopologicalDOF, fem::eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature>
			class FEMOperator {
			public:
				
				FEMOperator(const Assembler& assembler_, const mesh::Mesh& mesh_, const TopologicalDOF& topoDOF_, const Real time_, const Model& model_, const FormRegistry& forms_) : assembler(assembler_), mesh(mesh_), topoDOF(topoDOF_), time(time_), model(model_), forms(forms_) {}

				template<typename VectorType>
				void apply(const VectorType& x, VectorType& y) const {
					assembler.template assembleVector<TopologicalDOF::dofsPerNode, EvalEle, EvalQP, Model, FormRegistry, Quadrature>(mesh, topoDOF, time, model, forms, x, y);
				}

				Index size() const {
					return topoDOF.numFreeDOFs();
				}

			private:
				const Assembler& assembler;
				const mesh::Mesh& mesh;
				const TopologicalDOF& topoDOF;
				const Real time;
				const Model& model;
				const FormRegistry& forms;

			}; // class FEMOperator

		} // namespace op
	} // namespace linalg
} // namespace pdesolver

#endif
