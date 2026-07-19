#ifndef PDESOLVER_FEM_BOUNDARY_BOUNDARYCONDITION_HPP
#define PDESOLVER_FEM_BOUNDARY_BOUNDARYCONDITION_HPP

#include <concepts>

#include "core/Types.hpp"
#include "BoundaryCategory.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			template<typename Function>
			concept BoundaryFunction = requires(const Function f, Real time, const Real* x, Real* value) {
				{ Function::NumComponents } -> std::convertible_to<Index>;
				{ f.eval(time, x, value) } -> std::same_as<void>;
			}; // concept BoundaryFunction

			template<typename Function>
			struct BoundaryCondition {

				static constexpr Index NumComponents = Function::NumComponents;
				Int tag;
				BCCategory componentType[NumComponents];
				Function function;
			
			}; // struct BoundaryCondition

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
