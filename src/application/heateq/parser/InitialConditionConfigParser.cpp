#include "application/heateq/parser/InitialConditionConfigParser.hpp"

#include <stdexcept>

#include "solver/parser/NodalFieldReadConfigParser.hpp"
#include "solver/parser/NodalFieldWriteConfigParser.hpp"

pdesolver::application::heateq::config::InitialConditionConfig pdesolver::application::heateq::parser::InitialConditionConfigParser::parse(const YAML::Node& node) {

	pdesolver::application::heateq::config::InitialConditionConfig cfg;

	const YAML::Node& readNode = node["read"];
	if (!readNode) {
		throw std::runtime_error("InitialConditionConfigReader: missing required 'read' section in 'initial_condition'");
	}

	cfg.read = pdesolver::solver::parser::NodalFieldReadConfigParser::parse(readNode);
	cfg.write = pdesolver::solver::parser::NodalFieldWriteConfigParser::parse(node);

	return cfg;

}
