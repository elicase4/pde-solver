#include "application/heateq/BoundaryConditionConfigParser.hpp"
#include "application/heateq/ConductivityConfigParser.hpp"

#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::BoundaryConditionConfig::Type pdesolver::application::heateq::parser::BoundaryConditionConfigParser::parseBoundaryConditionType(const std::string& str) {

	if (str == "value") {
		return pdesolver::application::heateq::BoundaryConditionConfig::Type::Value;
	}

	if (str == "flux") {
		return pdesolver::application::heateq::BoundaryConditionConfig::Type::Flux;
	}

	throw std::runtime_error("Unknown boundary type: " + str);

}

pdesolver::application::heateq::config::BoundaryConditionConfig::Type pdesolver::application::heateq::parser::BoundaryConditionConfigParser::parseBoundaryConditionForm(const std::string& str) {

	if (str == "flux_bc") {
		return pdesolver::application::heateq::BoundaryConditionConfig::Form::FluxBC;
	}

	if (str == "value_bc") {
		return pdesolver::application::heateq::BoundaryConditionConfig::Form::ValueBC;
	}

	throw std::runtime_error("Unknown boundary form: " + str);

}

pdesolver::application::heateq::config::BoundaryConditionConfig pdesolver::application::heateq::parser::BoundaryConditionConfigParser::parse(const YAML::Node& node) {
	
	using io::YAMLReader;

	BoundaryConditionConfig cfg;

	cfg.boundaryID = YAMLReader::required<Int>(bc, "boundary");
	cfg.type = BoundaryConditionConfigParser::parseBoundaryConditionType(YAMLReader::required<std::string>(bc, "type"));
	cfg.expression = YAMLReader::required<std::string>(bc, "expression");
	for (auto& form : forms) {
		cfg.forms.push_back(BoundaryConditionConfigParser:parseBoundaryConditionForm(YAMLReader::required<std::string>(bc, "form")));
	}
	cfg.model = ConductivityConfigParser:parseConductivityType(YAMLReader::required<std::vector<std::string>>(bc, "model"));


	return cfg;

}
