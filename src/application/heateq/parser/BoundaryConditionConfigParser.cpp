#include "application/heateq/parser/BoundaryConditionConfigParser.hpp"
#include "application/heateq/parser/ConductivityConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::BoundaryConditionConfig::Type pdesolver::application::heateq::parser::BoundaryConditionConfigParser::parseBoundaryConditionType(const std::string& str) {

	if (str == "value") {
		return pdesolver::application::heateq::config::BoundaryConditionConfig::Type::Value;
	}

	if (str == "flux") {
		return pdesolver::application::heateq::config::BoundaryConditionConfig::Type::Flux;
	}

	throw std::runtime_error("Unknown boundary type: " + str);

}

pdesolver::application::heateq::config::BoundaryConditionConfig::Form pdesolver::application::heateq::parser::BoundaryConditionConfigParser::parseBoundaryConditionForm(const std::string& str) {

	if (str == "flux_bc") {
		return pdesolver::application::heateq::config::BoundaryConditionConfig::Form::FluxBC;
	}

	if (str == "value_bc") {
		return pdesolver::application::heateq::config::BoundaryConditionConfig::Form::ValueBC;
	}

	throw std::runtime_error("Unknown boundary form: " + str);

}

pdesolver::application::heateq::config::BoundaryConditionConfig pdesolver::application::heateq::parser::BoundaryConditionConfigParser::parse(const YAML::Node& node) {
	
	using io::YAMLReader;

	pdesolver::application::heateq::config::BoundaryConditionConfig cfg;

	cfg.boundaryID = YAMLReader::required<Int>(node, "boundary");
	cfg.type = BoundaryConditionConfigParser::parseBoundaryConditionType(YAMLReader::required<std::string>(node, "type"));

	// Value BCs take a single scalar expression; Flux BCs take one expression
	// per spatial component (a YAML sequence), since BoundaryFluxFunction
	// evaluates a full vector, not a scalar.
	if (cfg.type == pdesolver::application::heateq::config::BoundaryConditionConfig::Type::Flux) {
		cfg.fluxExpression = YAMLReader::required<std::vector<std::string>>(node, "expression");
	} else {
		cfg.expression = YAMLReader::required<std::string>(node, "expression");
	}
	for (auto& form : YAMLReader::required<std::vector<std::string>>(node, "forms")) {
		cfg.forms.push_back(BoundaryConditionConfigParser::parseBoundaryConditionForm(form));
	}
	cfg.model = pdesolver::application::heateq::parser::ConductivityConfigParser::parseConductivityType(YAMLReader::required<std::string>(node, "model"));


	return cfg;

}
