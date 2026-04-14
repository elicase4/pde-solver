#ifndef PDESOLVER_LINALG_CG_WORKSPACE_HPP
#define PDESOLVER_LINALG_CG_WORKSPACE_HPP

namespace pdesolver {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename VectorType>
					struct Workspace {

						VectorType r;
						VectorType p;
						VectorType Ap;
						VectorType z;

						explicit Workspace(Index n) : r(n), p(n), Ap(n), z(n) {}

					}; // struct Workspace

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace pdesolver

#endif
