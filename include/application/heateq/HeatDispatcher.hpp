#ifndef PDESOLVER_APPLICATION_HEATEQ_HEATDISPATCHER_HPP
#define PDESOLVER_APPLICATION_HEATEQ_HEATDISPATCHER_HPP

#include "application/heateq/config/HeatConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			class HeatDispatcher {
			public:

				static bool run(const config::HeatConfig& config);

			private:

				static bool dispatch2D(const config::HeatConfig& config);
				static bool dispatch3D(const config::HeatConfig& config);

			}; // class HeatDispatcher

		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
