#include "application/heateq/parser/HeatConfigParser.hpp"
#include "application/heateq/parser/BoundaryConditionConfigParser.hpp"
#include "application/heateq/parser/ConductivityConfigParser.hpp"
#include "application/heateq/parser/SourceConfigParser.hpp"

#include "solver/parser/DiscretizationConfigParser.hpp"
#include "solver/parser/MeshConfigParser.hpp"
#include "solver/parser/SolverConfigParser.hpp"
#include "solver/parser/OutputConfigParser.hpp"

#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::HeatConfig pdesolver::application::heateq::parser::HeatConfigParser::read(const std::string& filename) {

	using io::YAMLReader;

	const YAML::Node root = YAMLReader::loadFile(filename);

	pdesolver::application::heateq::config::HeatConfig cfg;

	const YAML::Node& mesh = root["mesh"];
	if (!mesh) {
		throw std::runtime_error("HeatConfigReader: missing required 'mesh' section in " + filename);
	}
	cfg.mesh = pdesolver::solver::parser::MeshConfigParser::parse(mesh);

	const YAML::Node& disc = root["discretization"];
	if (!disc) {
		throw std::runtime_error("HeatConfigReader: missing required 'discretization' section in " + filename);
	}
	cfg.discretization = pdesolver::solver::parser::DiscretizationConfigParser::parse(disc);

	const YAML::Node& solver = root["solver"];
	if (!solver) {
		throw std::runtime_error("HeatConfigParser: missing required 'solver' section in " + filename);
	}
	cfg.solver = pdesolver::solver::parser::SolverConfigParser::parse(solver);
	
	const YAML::Node& output = root["output"];
	if (!output) {
		throw std::runtime_error("HeatConfigParser: missing required 'output' section in " + filename);
	}
	cfg.output = pdesolver::solver::parser::OutputConfigParser::parse(output);
	
	const YAML::Node& bcs = root["boundary_conditions"];
	if (!bcs) {
		throw std::runtime_error("HeatConfigReader: missing required 'boundary_conditions' section in " + filename);
	}
	for (auto& bc : bcs) {
		cfg.boundaryConditions.push_back(pdesolver::application::heateq::parser::BoundaryConditionConfigParser::parse(bc));
	}

	const YAML::Node& physics = root["physics"];
	if (!physics) {
		throw std::runtime_error("HeatConfigReader: missing required 'physics' section in " + filename);
	}

	const YAML::Node& models = physics["models"];
	if (!models) {
		throw std::runtime_error("HeatConfigReader: missing required 'physics.models' section in " + filename);
	}

	const YAML::Node& cond = models["conductivity"];
	if (!cond) {
		throw std::runtime_error("HeatConfigReader: missing required 'physics.models.conductivity' section in " + filename);
	}
	cfg.conductivity = pdesolver::application::heateq::parser::ConductivityConfigParser::parse(cond);

	const YAML::Node& src = models["source"];
	if (!src) {
		throw std::runtime_error("HeatConfigReader: missing required 'physics.models.source' section in " + filename);
	}
	cfg.source = pdesolver::application::heateq::parser::SourceConfigParser::parse(src);
	
	const YAML::Node& backend = root["source"];
	if (!src) {
		throw std::runtime_error("HeatConfigReader: missing required 'physics.models.source' section in " + filename);
	}
	cfg.source = pdesolver::application::heateq::parser::SourceConfigParser::parse(src);

	return cfg;

}
