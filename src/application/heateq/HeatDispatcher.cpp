#include "application/heateq/HeatDispatcher.hpp"
#include "application/heateq/HeatStage.hpp"

#include "linalg/types/backend/CPU.hpp"

#include "fem/basis/LagrangeQuad.hpp"
#include "fem/basis/LagrangeHex.hpp"
#include "fem/quadrature/GaussQuadratureQuad.hpp"
#include "fem/quadrature/GaussQuadratureHex.hpp"
#include "fem/quadrature/GaussQuadrature1D.hpp"

#include "solver/driver/Steady.hpp"
#include "solver/driver/Transient.hpp"

#include <stdexcept>

namespace pdesolver {
	namespace application {
		namespace heateq {

			bool HeatDispatcher::runSteady(const HeatConfig& config) {

				switch (config.discretization.spatialDim) {
					case 2:  return dispatchSteady2D(config);
					case 3:  return dispatchSteady3D(config);
					default:
						throw std::runtime_error("HeatDispatcher: unsupported spatial_dim=" + std::to_string(config.discretization.spatialDim) + ". Valid: 2, 3");
				}
			}

			bool HeatDispatcher::dispatchSteady2D(const HeatConfig& config) {

				using Backend = linalg::types::backend::CPU;
				const Index Pxi  = config.discretization.basisOrderXi;
				const Index Peta = config.discretization.basisOrderEta;
				const Index Qxi  = config.discretization.quadraturePointXi;
				const Index Qeta = config.discretization.quadraturePointEta;

				if (Pxi != Peta) {
					throw std::runtime_error("HeatDispatcher: mixed basis orders (P_xi != P_eta) not yet supported");
				}

				// Dispatch on polynomial order.
				auto runP = [&]<Index P, Index Q>() -> bool {
					using Basis   = fem::basis::LagrangeQuad<P, P>;
					using QVol    = fem::quadrature::GaussQuadratureQuad<Q, Q>;
					using QBdy    = fem::quadrature::GaussQuadrature1D<Q>;
					using Stage   = HeatStage<Backend, Basis, QVol, QBdy>;
					using Driver  = solver::driver::Steady<Stage>;

					Stage  stage(config);
					Driver driver;
					return driver.solve(stage);
				};

				if (Pxi == 1 && Qxi == 2) return runP.template operator()<1, 2>();
				if (Pxi == 1 && Qxi == 3) return runP.template operator()<1, 3>();
				if (Pxi == 2 && Qxi == 3) return runP.template operator()<2, 3>();
				if (Pxi == 2 && Qxi == 4) return runP.template operator()<2, 4>();

				throw std::runtime_error("HeatDispatcher: unsupported 2D combination: P=" + std::to_string(Pxi) + " Q=" + std::to_string(Qxi) + ". Supported: (P=1,Q=2), (P=1,Q=3), (P=2,Q=3), (P=2,Q=4)");
			}

			bool HeatDispatcher::dispatchSteady3D(const HeatConfig& config) {

				using Backend = linalg::types::backend::CPU;
				const Index P = config.discretization.basisOrderXi;
				const Index Q = config.discretization.quadraturePointXi;

				if (P == 1 && Q == 2) {
					using Basis  = fem::basis::LagrangeHex<1, 1, 1>;
					using QVol   = fem::quadrature::GaussQuadratureHex<2, 2, 2>;
					using QBdy   = fem::quadrature::GaussQuadratureQuad<2, 2>;
					using Stage  = HeatStage<Backend, Basis, QVol, QBdy>;
					using Driver = solver::driver::Steady<Stage>;
					Stage  stage(config);
					Driver driver;
					return driver.solve(stage);
				}

				throw std::runtime_error("HeatDispatcher: unsupported 3D combination: P=" + std::to_string(P) + " Q=" + std::to_string(Q) + ". Supported: (P=1,Q=2)");

			}

			bool HeatDispatcher::runTransient(const HeatConfig& config) {

				switch (config.discretization.spatialDim) {
					case 2:  return dispatchTransient2D(config);
					case 3:  return dispatchTransient3D(config);
					default:
						throw std::runtime_error("HeatDispatcher: unsupported spatial_dim=" + std::to_string(config.discretization.spatialDim) );
				}
			}

			bool HeatDispatcher::dispatchTransient2D(const HeatConfig& config) {
				// TODO: wire up BackwardEuler / GeneralizedAlpha timestepper.
				// The pattern mirrors dispatchSteady2D:
				//   1. select Basis/QVol/QBdy template params from config
				//   2. construct a MassForm-aware HeatTransientStage (to be added)
				//   3. construct the TimeStepper from config.transient
				//   4. call Transient<Stage, Stepper, Vector>::solve(stage, stepper, U, U_prev)
				throw std::runtime_error("HeatDispatcher::dispatchTransient2D: not yet implemented");
			}

			bool HeatDispatcher::dispatchTransient3D(const HeatConfig& config) {
				throw std::runtime_error("HeatDispatcher::dispatchTransient3D: not yet implemented");
			}

		} // namespace heateq
	} // namespace application
} // namespace pdesolver
