#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>

namespace pde {
	namespace mesh{
		
		template<int DIM>
		class Mesh {
			public:
				
				// discretization parameters
				std::array<int, DIM> N;
				std::array<double, DIM> h;
				
				// domain bounds
				std::array<double, DIM> xmin;
				std::array<double, DIM> xmax;

				// constructors
				Mesh() = default;
				
				Mesh(const std::array<int, DIM>& N_, const std::array<double, DIM>& xmin_, const std::array<double, DIM>& xmax_): N(N_), xmin(xmin_), xmax(xmax_){
					
					// compute mesh spacing
					for (int d = 0; d < DIM; d++){
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

				std::size_t idx2flat(int i, int j, int k=0) const; 
				
				std::array<int,DIM> flat2idx(std::size_t flat_idx) const;

				std::array<double,DIM> idx2coord(int i, int j, int k=0) const;
				
				bool in_bounds(int i, int j, int k=0) const;

				// domain size 
				std::size_t size() const;
		};
	}
}
