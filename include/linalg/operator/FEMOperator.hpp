#ifndef PDESOLVER_LINALG_FEMOPERATOR_HPP
#define PDESOLVER_LINALG_FEMOPERATOR_HPP

#include "fem/assembly/Assembler.hpp"

namespace pdesolver {
	namespace linalg {
		namespace op {

			template<fem::eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form, typename Quadrature>
			class FEMOperator {
			public:
				
				FEMOperator(const fem::assembly::Assembler& assembler_, const mesh::Mesh& mesh_, const topology::topologicalDOF& topoDOF_, const Real time_, const Model& model_, const Form& form_) : assembler(assembler_), mesh(mesh_), topoDOF(topoDOF_), time(time_), model(model_), form(form_) {}

				template<typename VectorType>
				void apply(const VectorType& x, VectorType& y) const {
					assembler.assembleVector<EvalEle, EvalQP, Model, Form, Quadrature>(mesh, topoDOF, time, model, form, x, y);
				}

				Index size() const {
					return topoDOF.numFreeDofs();
				}

			private:
				fem::assembly::Assembler& assembler;
				mesh::Mesh& mesh;
				topology::topologicalDOF& topoDOF;
				Real time;
				Model model;
				Form form;

			}; // class FEMOperator

		} // namespace op
	} // namespace linalg
} // namespace pdesolver

#endif
