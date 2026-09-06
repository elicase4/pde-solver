#ifndef PDESOLVER_FEM_FORM_GATHERSELEMENTDATA_HPP
#define PDESOLVER_FEM_FORM_GATHERSELEMENTDATA_HPP

#include <concepts>

#include "core/Types.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {

			template<typename Form>
			concept GathersElementData = requires (const Form f, const Index* nodeIDs, Index nodesPerElement) {
				{ f.gatherElementData(nodeIDs, nodesPerElement) } -> std::same_as<void>;
			}; // concept GathersElementData

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
