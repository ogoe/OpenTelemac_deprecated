#CMake buildfile for sel2vtk
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
#
add_executable(sel2vtk sel2vtk.f90)
#
install(TARGETS sel2vtk
            RUNTIME DESTINATION bin)
#
#