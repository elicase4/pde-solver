#pragma once

#include "pde/grid.hpp"

namespace backend{
	
	template<int Dim>
	class Backend {
		
		public:
			virtual ~Backend() = default;

			virtual void applyOperator(const pde::Operator<Dim>& op, double* u, double* Au) = 0;
			
			virtual void updateSolution(const pde::Operator<Dim>& op, double* u, const double* rhs) = 0;
			
			virtual void computeResidual(const pde::Operator<Dim>& op, const double* u, const double* rhs, double* r) = 0;
	};

}
