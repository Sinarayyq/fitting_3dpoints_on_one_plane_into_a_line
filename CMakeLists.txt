# This CMakeList file is for project Surface_approximation.
# Please check if directories of libraries are correct before cmake.

project (iterative_ransac)

#configure_file (../../Packages/OpenNI2/Redist/OpenNI2.dll ../build COPYONLY)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 2.8.11)
if( POLICY CMP0043 )
  cmake_policy(SET CMP0043 OLD)
endif()

# Find PCL and add corrensponding include directories
find_package( PCL 1.8 REQUIRED )
include_directories( ${PCL_INCLUDE_DIRS} )
link_directories( ${PCL_LIBRARY_DIRS} )
add_definitions( ${PCL_DEFINITIONS} )


# Find OpenCV and add corrensponding include directories and additional libraries
set("OpenCV_DIR" ../../Packages/opencv/build)
find_package(OpenCV REQUIRED)
include_directories( ./include )


# Find CGAL and Qt5 and add corrensponding include directories
find_package(CGAL COMPONENTS Qt5)
include(${CGAL_USE_FILE})
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)
include_directories (BEFORE ../../include)
include_directories (BEFORE ./include)

# Include VCG library.
include_directories ( ../../Packages/vcglib-1.0.1 )

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND PCL_FOUND AND OpenCV_FOUND)

add_definitions(-DQT_NO_KEYWORDS)
add_definitions(-DNOMINMAX)

# UI files (Qt Designer files)
qt5_wrap_ui( DT_UI_FILES Alpha_shapes_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./Alpha_shapes_2.qrc )

# The executable itself.
add_executable ( iterative_ransac iterative_ransac.cpp)
qt5_use_modules (iterative_ransac Xml Script OpenGL Svg)
add_to_cached_list( CGAL_EXECUTABLE_TARGETS iterative_ransac )

# Link with Qt libraries
target_link_libraries (iterative_ransac ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries (iterative_ransac ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
# Link with PCL
target_link_libraries (iterative_ransac ${PCL_LIBRARIES})
# Link with OpenCV libraries
target_link_libraries (iterative_ransac ${OpenCV_LIBS})

else()

 
  message(STATUS "NOTICE: iterative_ransac project requires PCL, OpenCV, CGAL and Qt5, and will not be compiled.")

endif()

