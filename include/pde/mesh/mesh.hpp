#pragma once
#include <array>
#include <cstddef>
#include <stdexcept>

namespace pde {
	namespace mesh{
		
		class Mesh {
			public:
				
				// discretization parameters
				std::array<int, 3> N;
				std::array<double, 3> h;
				
				// domain bounds
				std::array<double, 3> xmin;
				std::array<double, 3> xmax;

				// constructors
				Mesh() = default;
				
				Mesh(const std::array<int, 3>& N_, const std::array<int, 3>& xmin_, const std::array<int, 3>& xmax_): N(N_), xmin(xmin_), xmax(xmax_){
					
					// compute mesh spacing
					for (int d = 0; d < 3; d++){
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
				inline std::size_t idx(int i, int j, int k) const {
					return (std::size_t(i*N[1]) + std::size_t(j) + std::size_t(k*N[0]*N[1]));
				}

				// check bounds
				inline bool in_bounds(int i, int j, int k) const {
					return ((i >= 0 && i < N[0]) && (j >= 0 && j < N[1]) && (k >= 0 && k < N[2]));
				}
				
				// indexing to coordinate
				inline std::array<double,3> coord(int i, int j, int k) const {
					std::array<double,3> out {};
					out[0] = xmin[0] + (i + 1)*h[0];
					out[1] = xmin[1] + (j + 1)*h[1];
					out[2] = xmin[2] + (k + 1)*h[2];
					return out;
				}

				// domain size 
				std::size_t size() const {
					return std::size_t(N[0] * N[1] * N[2]);
				}
		};
	}
}
