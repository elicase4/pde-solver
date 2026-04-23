#ifndef PDESOLVER_SOLVER_PRECONDITIONER_IDENTITY_HPP
#define PDESOLVER_SOLVER_PRECONDITIONER_IDENTITY_HPP

#include "linalg/operations/VectorOps.hpp"

namespace pdesolver {
	namespace linalg {
		namespace solver {
			namespace preconditioner {

				template<typename VectorType>
				class Identity {
				public:
					
					void apply(const VectorType& r, VectorType& z){
						operations::copy(r, z); // z = r
					}

				}; // class Identity

			} // namespaace preconditioner
		} // namespaace solver
	} // namespaace linalg
} // namespaace pdesolver

#endif
