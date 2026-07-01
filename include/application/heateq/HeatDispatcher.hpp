#ifndef PDESOLVER_APPLICATION_HEATEQ_HEATDISPATCHER_HPP
#define PDESOLVER_APPLICATION_HEATEQ_HEATDISPATCHER_HPP

#include "application/heateq/HeatConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			class HeatDispatcher {
			public:

				static bool runSteady(const HeatConfig& config);

				static bool runTransient(const HeatConfig& config);

			private:

				static bool dispatchSteady2D(const HeatConfig& config);
				static bool dispatchSteady3D(const HeatConfig& config);

				static bool dispatchTransient2D(const HeatConfig& config);
				static bool dispatchTransient3D(const HeatConfig& config);

			}; // class HeatDispatcher

		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
