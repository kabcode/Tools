cmake_minimum_required(VERSION 3.10)
project(InterpolationComparison)

set(CMAKE_CXX_STANDARD 17)
find_package(ITK REQUIRED
	COMPONENTS
		RTK
		ITKCommon
		ITKLabelMap
		ITKImageIntensity
		ITKThresholding
		ITKIOImageBase
		ITKIOMeta
		ITKIONRRD
)

include(${ITK_USE_FILE})

set (SRCS
	${CMAKE_PROJECT_NAME}.cxx
	HighPrecisionTimer.h
)

add_executable(${CMAKE_PROJECT_NAME}
	${SRCS}
)

target_link_libraries(${CMAKE_PROJECT_NAME}
	${ITK_LIBRARIES}	
)