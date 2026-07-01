cmake_minimum_required(VERSION 3.18)

find_package(exprtk REQUIRED)

add_executable(expression
	expression/ScalarExpression.cpp
	expression/VectorExpression.cpp
	expression/tensorExpression.cpp
)
