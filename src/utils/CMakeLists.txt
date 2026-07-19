find_path(EXPRTK_INCLUDE_DIR
    NAMES exprtk.hpp
    PATH_SUFFIXES include
)

if (NOT EXPRTK_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find exprtk.hpp")
endif()

add_library(pde_expression
	expression/ScalarExpression.cpp
	expression/VectorExpression.cpp
	expression/TensorExpression.cpp
)

target_include_directories(pde_expression
    PUBLIC
        ${EXPRTK_INCLUDE_DIR}
)
