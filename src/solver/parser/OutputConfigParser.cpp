#include "solver/config/OutputConfigParser.hpp"
#include "io/YAML/Reader.hpp"

pdesolver::solver::config::OutputConfig pdesolver::solver::parser::OutputConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	OutputConfig cfg;

	cfg.directory = YAMLReader::optional<std::string>(node, "directory", "output");
	cfg.vtk = YAMLReader::optional<bool>(node, "vtk", true);
	cfg.writeFrequency = YAMLReader::optional<Index>(node, "write_frequency", 1);
	cfg.prefix = YAMLReader::optional<std::string>(node, "prefix", "solution");

	return cfg;

}
