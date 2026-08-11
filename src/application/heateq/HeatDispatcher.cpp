#include <iostream>
#include <stdexcept>
#include <utility>

#include "application/heateq/HeatDispatcher.hpp"
#include "application/heateq/problem/HeatProblem.hpp"
#include "application/heateq/stage/HeatStage.hpp"

#include "equations/heateq/HeatEquation.hpp"

#include "fem/dispatch/DiscretizationDispatch.hpp"

#include "io/MeshIO.hpp"

#include "linalg/types/backend/CPU.hpp"

#include "mesh/Mesh.hpp"

#include "solver/driver/Steady.hpp"

bool pdesolver::application::heateq::HeatDispatcher::run(const pdesolver::application::heateq::config::HeatConfig& config) {

	namespace sconfig = pdesolver::solver::config;

	if (config.backend.type != sconfig::BackendConfig::Type::CPU) {
		throw std::runtime_error("HeatDispatcher: only the CPU backend is supported so far");
	}

	if (config.solver.driver.type != sconfig::DriverConfig::Type::Steady) {
		throw std::runtime_error("HeatDispatcher: only the steady driver is supported so far");
	}

	if (!config.solver.linear.has_value()) {
		throw std::runtime_error("HeatDispatcher: solver.linear config is required");
	}

	pdesolver::mesh::Mesh mesh;
	pdesolver::io::MeshIO::readBinary(mesh, config.mesh.file);

	if (mesh.data.basisType != pdesolver::mesh::BasisType::Lagrange) {
		throw std::runtime_error("HeatDispatcher: only the Lagrange basis is supported so far");
	}

	const Index nsd = mesh.data.spatialDim;
	const Index npd = mesh.data.parametricDim;
	const auto family = mesh.data.elementFamily;

	const Index basisOrderX = mesh.data.basisOrder.size() > 0 ? mesh.data.basisOrder[0] : 1;
	const Index basisOrderY = mesh.data.basisOrder.size() > 1 ? mesh.data.basisOrder[1] : 1;
	const Index basisOrderZ = mesh.data.basisOrder.size() > 2 ? mesh.data.basisOrder[2] : 1;

	std::cout << "[heateq] mesh: " << mesh.data.numElements << " elements, nsd=" << nsd << " npd=" << npd << "\n";

	bool converged = false;

	auto visitor = [&]<Index NSD, typename BasisT, typename QuadVolT, typename QuadBdyT>() -> bool {

		using HeatEqBundle = pdesolver::equations::HeatEquation<NSD, BasisT, QuadVolT, QuadBdyT>;
		using ProblemT = pdesolver::application::heateq::problem::HeatProblem<pdesolver::linalg::types::backend::CPU, HeatEqBundle>;
		using StageT = pdesolver::application::heateq::stage::HeatStage<ProblemT>;

		ProblemT heatProblem(config, std::move(mesh));
		StageT stage(heatProblem);

		pdesolver::solver::driver::Steady<StageT> driver;
		converged = driver.solve(stage);

		heatProblem.writeOutput(0);
		heatProblem.writeLog();

		std::cout << "[heateq] " << (converged ? "converged" : "did not converge") << "\n";

		return converged;

	};

	const auto& d = config.discretization;
	const bool matched = pdesolver::fem::dispatch::dispatch(nsd, npd, family, basisOrderX, basisOrderY, basisOrderZ, d.quadrature.xi, d.quadrature.eta, d.quadrature.zeta, visitor);

	if (!matched) {
		throw std::runtime_error("HeatDispatcher: unsupported (nsd, npd, basis order, quadrature order) combination for this mesh/config");
	}

	return converged;

}
