#ifndef PDESOLVER_FEM_QUANTITY_QUANTITYUNITS_HPP
#define PDESOLVER_FEM_QUANTITY_QUANTITYUNITS_HPP

#include <string>

#include "core/Types.hpp"

#include "fem/quantity/Reduction.hpp"

namespace pdesolver::fem::quantity {

	template<typename Form>
	std::string unitFor(Reduction mode, Index boundaryDim) {

		std::string unit = Form::BaseUnit;
		if (mode == Reduction::Average) {
			unit += "/m^" + std::to_string(boundaryDim);
		}

		return unit;

	}

} // namespace pdesolver::fem::quantity

#endif
