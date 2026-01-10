#ifndef PDESOLVER_BOUNDARYCONDITION_HPP
#define PDESOLVER_BOUNDARYCONDITION_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			enum class BCType {
				Essential,
				Natural
			}; // class  BCType
		

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver
