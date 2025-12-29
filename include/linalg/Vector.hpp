#ifndef PDESOLVER_VECTOR_HPP
#define PDESOLVER_VECTOR_HPP

#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace linalg {
		
		class Vector {
		public:
			Vector(Index size = 0) : data(size, 0.0) {}

			// Access
			Real& operator[](Index i) { return data[i]; }
			Real operator[](Index i) const { return data[i]; }

			// Size
			Index size() const { return data.size(); }
			void resize(Index n) const { return data.resize(n, 0.0); }

			// Data Access
			Real* getData() { return data.data(); }
			const Real* getData() const { return data.data(); }

			// Utilities
			void setZero() { for (auto& x : data) x = 0.0; }
			
		private:
			std::vector<Real> data;
		
		}; // class Vector
		
	} // namespace linalg
} // namespace pdesolver

#endif
