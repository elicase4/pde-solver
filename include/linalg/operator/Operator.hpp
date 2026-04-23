#ifndef PDESOLVER_LINALG_OPERATOR_HPP
#define PDESOLVER_LINALG_OPERATOR_HPP

#include <concepts>

namespace pdesolver {
	namespace linalg {
		namespace op {

			template<typename Operator, typename VectorType>
			concept LinearOperator = requires(const Operator A, const VectorType& x, VectorType& y){
				{ A.apply(x, y) } -> std::same_as<void>;
			};

		} // namespace op
	} // namespace linalg
} // namespace pdesolver

#endif
