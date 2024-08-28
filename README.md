# SingletEWPT

Program that does stuff

### Requirements

Requires C++17 compliant compiler and CMake 3.15 or newer. We depend on GNU Scientific Library (GSL, https://www.gnu.org/software/gsl/) and ```dlib``` (http://dlib.net). You will need to manually install the GSL libraries. For ```dlib``` you can either manually install the libraries, or let CMake automatically fetch it for you during build stage (see below).

### Compiling and installing

```
cmake -B build
cmake --build build
cmake --install build
```
The binaries get installed to ```bin/``` in the project root directory. 

For a debug build, configure with the ```-DCMAKE_BUILD_TYPE=Debug``` flag.
