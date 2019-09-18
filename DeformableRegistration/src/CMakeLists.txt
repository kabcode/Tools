cmake_minimum_required(VERSION 3.11)
project(DeformableRegistration)

find_package(ITK REQUIRED COMPONENTS ITKImageCompare ITKRegistrationMethodsv4 ITKIOImageBase ITKIOTransformBase ITKIOXML ITKImageGrid ITKImageIntensity ITKThresholding ITKTransform ITKImageIO ITKTransformIO)
include(${ITK_USE_FILE})
set(SRC
	DeformableRegistration.cxx
)
add_executable(${CMAKE_PROJECT_NAME} ${SRC})
target_link_libraries( ${CMAKE_PROJECT_NAME}
	${ITK_LIBRARIES}
)
