### OpenTelemac
This is the [Telemac modelling system](http://www.opentelemac.org/ "OpenTelemac") with all bugfixes from the official svn repository and some enhancements made by [Oliver Goethel](http://github.com/ogoe "ogoe").

So far the main changes compared to the official version are:

* CMake build system
* changed directory structure for easier diff/merge/branch operations
* adapted perl scripts and dicos for changed directory structure
* additional sediment transport equation(s) 
* alternative sandslide algorithm (only in the 6.0 branch)
* wave module for Telemac 3D (Airy, Stokes and Stream Function) (not in the repo yet)
* ....


If you like to have the original (svn) version of the Telemac system, you might want to take a look at the [svn mirror repository](http://github.com/ogoe/OpenTelemac-svn-mirror "Telemac svn Mirror").

### Build instructions

This (repository) version comes with CMake build files. I strongly recommend to do an out of source build. So the only thing to do (assuming you have already installed [cmake](http://www.cmake.org)) to build from source is:

> `git clone http://github.com/ogoe/OpenTelemac.git opentelemac`
>
> `cd opentelemac`
> 
>`mkdir build && cd build`
>
>`cmake ..`
>
>`make && make install`

CMake will recognize if you have a (Fortran) MPI installation on your system and build the parallel library and executables accordingly. For `partel` to be build, the `libmetis.a` library has to be available in your `parallel` directory.

Enjoy.

