cmake_minimum_required(VERSION 3.10)
project(TransformVolumeToDeviceCOS)

find_package(ITK REQUIRED
	COMPONENTS
		ITKIOImageBase
		ITKIOTransformBase
		ITKIOXML
		ITKImageGrid
		ITKImageIntensity
		ITKThresholding
		ITKTransform
		ITKImageIO
		ITKTransformIO
)

include(${ITK_USE_FILE})
set(SRC
	${CMAKE_PROJECT_NAME}.cxx
)
add_executable(${CMAKE_PROJECT_NAME} ${SRC})
target_link_libraries( ${CMAKE_PROJECT_NAME}
	${ITK_LIBRARIES}
)