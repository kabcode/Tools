cmake_minimum_required(VERSION 3.10)
project(OrganMotion)

set(CMAKE_CXX_STANDARD 17)

find_package(ITK REQUIRED COMPONENTS ITKIOImageBase ITKIOTransformBase ITKImageGrid ITKImageIntensity ITKTransform ITKImageIO)
include(${ITK_USE_FILE})
set(SRC
	OrganMotion.cxx
)
add_executable(${CMAKE_PROJECT_NAME} ${SRC})
target_link_libraries( ${CMAKE_PROJECT_NAME}
	${ITK_LIBRARIES}
)