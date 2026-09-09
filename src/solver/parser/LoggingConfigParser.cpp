#include "solver/parser/LoggingConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"

pdesolver::solver::config::LoggerConfig::Type pdesolver::solver::parser::LoggingConfigParser::parseLoggerType(const std::string& str) {

	if (str == "console") {
		return pdesolver::solver::config::LoggerConfig::Type::Console;
	}

	if (str == "none") {
		return pdesolver::solver::config::LoggerConfig::Type::None;
	}

	throw std::runtime_error("Unknown logger type: " + str);

}

pdesolver::solver::config::LoggingConfig pdesolver::solver::parser::LoggingConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::solver::config::LoggingConfig cfg;

	if (!node) {
		return cfg;
	}

	const YAML::Node& solverNode = node["solver"];
	if (solverNode) {
		cfg.solver.type = LoggingConfigParser::parseLoggerType(YAMLReader::optional<std::string>(solverNode, "type", "console"));
		cfg.solver.textFile = YAMLReader::optional<std::string>(solverNode, "text_file", "");
		cfg.solver.csvFile = YAMLReader::optional<std::string>(solverNode, "csv_file", "");
	}

	const YAML::Node& driverNode = node["driver"];
	if (driverNode) {
		cfg.driver.type = LoggingConfigParser::parseLoggerType(YAMLReader::optional<std::string>(driverNode, "type", "console"));
		cfg.driver.textFile = YAMLReader::optional<std::string>(driverNode, "text_file", "");
	}

	return cfg;

}
