#ifndef PDESOLVER_LINEARASSEMBLER_HPP
#define PDESOLVER_LINEARASSEMBLER_HPP

#include "core/Types.hpp"
#include "mesh/MeshBase.hpp"
#include "linalg/SparseMatrix.hpp"
#include "linalg/Vector.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {
			
			template<typename Element, typename Quadrature, typename Backend>
			class LinearAssembler{
			public:
				
				// allocation functions
				static linalg::SparseMatrix createMatrix(const mesh::MeshBase& mesh);
				static linalg::Vector createVector(const mesh::MeshBase& mesh);
				
				// matrix assembly
				template<typename BilinearForm>
				static void assmebleMatrix(const mesh::MeshBase& mesh, linalg::SparseMatrix& K);
				
				// operator assembly
				template<typename BilinearForm>
				static void assmebleOperator(const mesh::MeshBase& mesh, linalg::Vector& o);
				
				// vector assembly
				template<typename LinearForm>
				static void assmebleMatrix(const mesh::MeshBase& mesh, linalg::Vector& f);

				// deallocation functions
				static linalg::SparseMatrix deleteMatrix();
				static linalg::Vector deleteVector();
			
			}; // class LinearAssember

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#endif
