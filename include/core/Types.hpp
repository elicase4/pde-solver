#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstddef>
#include <vector>
#include <array>

namespace pdesolver {
	
	// floating point types
	using Real = double;

	// integer types
	using Index = std::size_t;
	using Int = int;

	// array types
	template<typename T, Int N>
	using Array = std::array<T, N>;
	
	using Point2D = Array<Real, 2>;
	using Point3D = Array<Real, 3>;
}

#endif
