#pragma once
#include <array>
#include <cstddef>
#include <stdexcept>


namespace pde {
	namespace mesh{
		
		const int DIM = 2;
		
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

				// cartesian to flat indexing
				inline std::size_t idx(int i, int j) const {
					return (std::size_t(i*N[1]) + std::size_t(j));
				}

				// check bounds
				inline bool in_bounds(int i, int j) const {
					return ((i >= 0 && i < N[0]) && (j >= 0 && j < N[1]));
				}
				
				// indexing to coordinate
				inline std::array<double,DIM> coord(int i, int j) const {
					std::array<double,DIM> out {};
					out[0] = xmin[0] + (i + 1)*h[0];
					out[1] = xmin[1] + (j + 1)*h[1];
					return out;
				}

				// domain size 
				std::size_t size() const {
					return std::size_t(N[0] * N[1]);
				}
		};
	}
}
