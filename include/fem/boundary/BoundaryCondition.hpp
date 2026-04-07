#ifndef PDESOLVER_BOUNDARYCONDITION_HPP
#define PDESOLVER_BOUNDARYCONDITION_HPP

#include <utility>

#include "core/Types.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			enum class BCCategory {
				None,
				Essential,
				Natural
			}; // enum class BCCategory

			template<typename Function>
			concept BoundaryFunction = requires (const Function f, const Real time, const Real* x, Real* outValue) {
				
				{ Function::NumComponents } -> std::convertible_to<Index>;
				{ f.eval(time, x, outValue) } -> std::same_as<void>;

			}; // concept BoundaryFunction

			template<typename Function>
			struct BoundaryCondition {
				
				static constexpr Index NumComponents = Function::NumComponents;
				
				Int tag;
				BCCategory componentType[Function::NumComponents];
				Function f;

			}; // struct BoundaryCondition

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
