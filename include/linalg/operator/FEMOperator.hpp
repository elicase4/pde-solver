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

			template<typename Assembler, fem::eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form, typename Quadrature>
			class FEMOperator {
			public:
				
				FEMOperator(const Assembler& assembler_, const mesh::Mesh& mesh_, const topology::TopologicalDOF& topoDOF_, const Real time_, const Model& model_, const Form& form_) : assembler(assembler_), mesh(mesh_), topoDOF(topoDOF_), time(time_), model(model_), form(form_) {}

				template<typename VectorType>
				void apply(const VectorType& x, VectorType& y) const {
					assembler.template assembleVector<EvalEle, EvalQP, Model, Form, Quadrature>(mesh, topoDOF, time, model, form, x, y);
				}

				Index size() const {
					return topoDOF.numFreeDOFs();
				}

			private:
				const Assembler& assembler;
				const mesh::Mesh& mesh;
				const topology::TopologicalDOF& topoDOF;
				const Real time;
				const Model model;
				const Form form;

			}; // class FEMOperator

		} // namespace op
	} // namespace linalg
} // namespace pdesolver

#endif
