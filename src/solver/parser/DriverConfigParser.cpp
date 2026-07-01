#include "solver/config/DriverConfigParser.hpp"
#include "io/YAML/Reader.hpp"

pdesolver::solver::config::DriverConfig::Type pdesolver::solver::parser::DriverConfigParser::parseDriverType(const std::string& str) {

	if (str == "steady")
		return DriverConfig::Type::Steady;
	if (str == "transient")
		return DriverConfig::Type::Transient;

	throw std::runtime_error("Unknown driver type: '" + str + "'. Valid options: steady, transient");

}

pdesolver::solver::config::DriverConfig pdesolver::solver::parser::DriverConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	SolverConfig cfg;

	cfg.driver = parseDriverType(YAMLReader::required<std::string>(node, "driver"));

	return cfg;

}
