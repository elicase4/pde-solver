#ifndef PDESOLVER_APPLICATION_HEATEQ_HEATAPPLICATION_HPP
#define PDESOLVER_APPLICATION_HEATEQ_HEATAPPLICATION_HPP

#include "application/heateq/config/HeatConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			class HeatApplication {
			public:

				explicit HeatApplication(const config::HeatConfig& config);

				int run();

			private:

				config::HeatConfig config_;

			}; // class HeatApplication

		}
	}
}

#endif
