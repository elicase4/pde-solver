#pragma once
#include <array>
#include <cstddef>
#include <stdexcept>

namespace pde {
	
	template<int Dim>
	class Grid {
		public:
			
			static_assert(Dim==2 || Dim==3, "Input Dimension must be 2 or 3");
			
			std::array<int, Dim> N; // interior points per axis
			std::array<double, Dim> h; // node spacing per axis
			
			// domain bounds
			std::array<double,Dim> xmin;
			std::array<double,Dim> xmax;
			
			// constructors
			Grid() = default;
			
			Grid(const std::array<int, Dim>& N_, const std::array<int, Dim>& xmin_, const std::array<int, Dim>& xmax_): N(N_), xmin(xmin_), xmax(xmax_){
				
				for (int d = 0; d < Dim; d++){
					if (N[d] <= 0){
						throw std::runtime_error("Grid dimension N[d] must be positive.");
					}
					if (xmax[d] <= xmin[d]){
						throw std::runtime_error("Grid domain bounds are invalid.");
					}
					// store interior spacing
					h[d] = (xmax[d] - xmin[d]) / (N[d] + 1);
				}
			
			}
			
			inline std::size_t idx(int i, int j, int k=0) const; // cartesian to flat indexing
			inline bool in_bounds(int i, int j, int k=0) const; // check in bounds indexing
			inline std::array<double,Dim> coord(int i, int j, int k=0) const; // get index to coordinates

			// domain size 
			std::size_t size() const {
				std::size_t s = 1;
				for (int d=0; d<Dim; d++){
					s *= static_cast<std::size_t>(N[d]);
				}
				return s;
			}
	};
}
