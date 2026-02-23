#ifndef PDESOLVER_BOUNDARYCONDITION_HPP
#define PDESOLVER_BOUNDARYCONDITION_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			enum class BCCategory {
				Essential,
				Natural
			}; // enum class BCCategory

			template<typename BC>
			concept BoundaryCondition = requires (const Real t, const Real* x, Real* v) {
				{ BC::evaluate(t, x, v) } -> std::same_as<void>;
				{ BC::category() } -> std::same_as<BCCategory>;
				{ BC::numComponents() } -> std::convertible_to<Index>;
			}; // concept BoundaryCondition

			template<typename BC>
			concept EssentialBC = BoundaryCondition<BC> && (BC::category() == BCCategory::Essential);

			template<typename BC>
			concept NaturalBC = BoundaryCondition<BC> && (BC::category() == BCCategory::Natural);

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
