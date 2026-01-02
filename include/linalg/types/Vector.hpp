#ifndef PDESOLVER_VECTOR_HPP
#define PDESOLVER_VECTOR_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace linalg {
		namespace types {

			template<typename T, typename Backend>
			class Vector {
			public:
				explicit Vector(Index size) : size_(size), data_(Backend::template alloc<T>(size)) {};
				
				// Move only operations
				Vector(const Vector&) = delete;
				Vector& operator=(const Vector&) = delete;
				Vector(Vector&&) = default;
				Vector& operator=(Vector&&) = default;

				// Size
				Index size() const { return size_; }

				// Access
				T* data() { return data_.get(); }
				const T* data() const { return data_.get(); }
			
			private:
				Index size_;
				typename Backend::template Ptr<T> data_;

			}; // class Vector
		
			} // namespace types
		} // namespace linalg
	} // namespace pdesolver

#endif
