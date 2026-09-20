#ifndef PDESOLVER_FEM_FORM_SCALEDFORM_HPP
#define PDESOLVER_FEM_FORM_SCALEDFORM_HPP

#include <cstring>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"

namespace pdesolver {
	namespace fem {
		namespace form {

			template<Index numDOFs, typename Form>
			class ScaledForm {
			public:

				constexpr ScaledForm() = default;
				constexpr explicit ScaledForm(Real coefficient) : coefficient_(coefficient) {}

				template<typename QuadraturePoint>
				PDE_HOST PDE_DEVICE void computeElementLevelMatrix(const QuadraturePoint& qp, Real* Ke) const {

					constexpr Index maxElementDOFs = fem::dispatch::kMaxNodesPerElement<QuadraturePoint::ParametricDim> * numDOFs;
					Real buffer[maxElementDOFs * maxElementDOFs];
					std::memset(buffer, 0.0, sizeof(buffer));

					inner_.computeElementLevelMatrix(qp, buffer);

					const Index n = qp.nodesPerElement() * numDOFs;
					for (Index i = 0; i < n * n; ++i) {
						Ke[i] += coefficient_ * buffer[i];
					}

				}

				template<typename QuadraturePoint>
				PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePoint& qp, Real* Ue, Real* Oe) const {

					Real buffer[fem::dispatch::kMaxNodesPerElement<QuadraturePoint::ParametricDim> * numDOFs];
					std::memset(buffer, 0.0, sizeof(buffer));

					inner_.computeElementLevelVector(qp, Ue, buffer);

					const Index n = qp.nodesPerElement() * numDOFs;
					for (Index i = 0; i < n; ++i) {
						Oe[i] += coefficient_ * buffer[i];
					}

				}

			private:

				Form inner_;
				Real coefficient_ = Real(1);

			}; // class ScaledForm

		} // namespace form
	} // namespace fem
} // namespace pdesolver

#endif
