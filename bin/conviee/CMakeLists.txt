#CMake buildfile for conviee
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
#
add_executable(conviee.exe conviee.f)
#
install(TARGETS conviee.exe
            RUNTIME DESTINATION bin)
#
#
