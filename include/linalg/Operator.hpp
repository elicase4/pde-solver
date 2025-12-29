#ifndef PDESOLVER_OPERATOR_HPP
#define PDESOLVER_OPERATOR_HPP

#include "linalg/Vector.hpp"

namespace pdesolver {
	namespace linalg {
		
		class Operator {
		public:
			Operator(Index nRows_) : data(nRows, 0.0), nRows(nRows_) {};
		private:
			Index nRows;
			std::vector<Real> data;
		}; // class Operator

	} // namespace linalg
} // namespace pdesolver

#endif
