#ifndef PDESOLVER_ASSEMBLER_HPP
#define PDESOLVER_ASSEMBLER_HPP

#include <algorithm>

#include "config/Platform.hpp"

#include "core/Types.hpp"

#include "fem/eval/EvalElement.hpp"
#include "fem/eval/EvalQuadraturePoint.hpp"
#include "fem/eval/EvalModel.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/form/LinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"
#include "fem/form/NonlinearForm.hpp"

#include "mesh/Mesh.hpp"

#include "linalg/types/Matrix.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/types/DistributedCSRMatrix.hpp"
#include "linalg/types/DistributedVector.hpp"

#include "topology/TopologicalDOF.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {

			namespace eval = pdesolver::fem::eval;

			template<typename Backend>
			class Assembler {
			public:
			
				// constructor
				Assembler(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF) : mesh_(mesh), topoDOF_(topoDOF) {}

				// allocation function
				linalg::types::CSRMatrix<Real, Backend> createMatrix();
				
				// allocation function
				linalg::types::Vector<Real, Backend> createVector();
				
				// matrix assembly
				template<eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form, typename Quadrature>
				requires eval::EvalQuadraturePoint<EvalQP, EvalEle> && eval::EvalModel<Model, EvalQP>
				void assembleMatrix(const Real time, const Model& model, const Form& form, const linalg::types::Vector<Real, Backend>& U, linalg::types::CSRMatrix<Real, Backend>& K);
				
				// vector assembly
				template<eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form, typename Quadrature>
				requires eval::EvalQuadraturePoint<EvalQP, EvalEle> && eval::EvalModel<Model, EvalQP>
				void assembleVector(const Real time, const Model& model, const Form& form, const linalg::types::Vector<Real, Backend>& U, linalg::types::Vector<Real, Backend>& F);

			private:
				const mesh::Mesh& mesh_;
				const topology::TopologicalDOF& topoDOF_;

			}; // class Assembler

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/Assembler.tpp"

#endif
