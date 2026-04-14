#ifndef PDESOLVER_LINALG_CG_CONFIG_HPP
#define PDESOLVER_LINALG_CG_CONFIG_HPP

namespace pdesolver {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename VectorType>
					struct Config {
						using DataType = typename VectorType::value_type;

						DataType tol = 1e-8;
						Index maxIters = 1000;
					}; // struct Config

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace pdesolver

#endif
