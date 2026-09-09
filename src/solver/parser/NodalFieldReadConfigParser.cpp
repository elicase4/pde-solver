#include "solver/parser/NodalFieldReadConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"

pdesolver::solver::config::NodalFieldReadConfig::Mode pdesolver::solver::parser::NodalFieldReadConfigParser::parseMode(const std::string& str) {

	if (str == "expression") {
		return pdesolver::solver::config::NodalFieldReadConfig::Mode::Expression;
	}

	if (str == "file") {
		return pdesolver::solver::config::NodalFieldReadConfig::Mode::File;
	}

	throw std::runtime_error("Unknown read mode: " + str);

}

pdesolver::solver::config::NodalFieldReadConfig pdesolver::solver::parser::NodalFieldReadConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::solver::config::NodalFieldReadConfig cfg;

	cfg.mode = NodalFieldReadConfigParser::parseMode(YAMLReader::optional<std::string>(node, "mode", "expression"));

	if (cfg.mode == pdesolver::solver::config::NodalFieldReadConfig::Mode::File) {
		cfg.file = YAMLReader::required<std::string>(node, "file");
	} else {
		cfg.expression = YAMLReader::required<std::string>(node, "expression");
	}

	cfg.unit = YAMLReader::required<std::string>(node, "unit");

	return cfg;

}
