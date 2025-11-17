#pragma once

#include <array>
#include <functional>
#include <stdexcept>

#include "pde/mesh/mesh.hpp"
#include "pde/solution/solution.hpp"
#include "pde/rhs/rhs.hpp"

namespace pde {
	namespace boundary {
		
		template<int DIM>
		class Dirichlet {
			public:
				Dirichlet(const pde::mesh::Mesh<DIM>& mesh_, const pde::solution::Solution<DIM>& sol_, const pde::rhs::RHS<DIM>& rhs_) : mesh(mesh_), solution(sol_), rhs(rhs_) {};

				void set(pde::mesh::BoundaryFace face, std::function<double(const std::array<double, DIM>&)> f){
					bc_funcs[(int) face] = f;
				}

				void eval(double* u, double* f) const{
					
					// add contribution to u
					std::ptrdiff_t sol_flat = u - solution.buffer_base;
					std::array<int,DIM> sol_idx = mesh.flat2idx(static_cast<std::size_t>(sol_flat));
					
					pde::mesh::BoundaryFace face = mesh.get_boundary_face(sol_idx[0], sol_idx[1], sol_idx[2]);
					if(face == pde::mesh::BoundaryFace::None){
						throw std::runtime_error("Boundary function evaluation used on interior node.");
					}
					
					auto sol_coord = mesh.idx2coord(sol_idx[0], sol_idx[1], sol_idx[2]);
					
					u[0] = bc_funcs[(int) face](sol_coord);
					
					// add contribution to rhs
					int axis = ((int) face) / 2;
					int offset = ((int) face % 2 == 0) ? 1 : -1;
					
					std::array<int,DIM> rhs_idx = {0};
					rhs_idx[axis] += offset;
					std::size_t rhs_flat = mesh.idx2flat(rhs_idx[0], rhs_idx[1], rhs_idx[2]);
					
					f[rhs_flat] -= (bc_funcs[(int) face](sol_coord)) / (mesh.h[axis] * mesh.h[axis]);
				}

			private:
				const pde::mesh::Mesh<DIM>& mesh;
				const pde::solution::Solution<DIM>& solution;
				const pde::rhs::RHS<DIM>& rhs;
				std::array<std::function<double(const std::array<double,DIM>&)>, 2*DIM> bc_funcs{};
		};
	}
}
