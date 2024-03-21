# SingletEWPT

Program that does stuff

## Dependencies

- dlib C++ library
- GSL

## Installation

### Installing dependencies

You will need to manually install the GSL libraries.

For dlib, you can either manually install the libraries from http://dlib.net, or let CMake automatically fetch it for you in the next step.


### Building the program

```
cmake -B build
cmake --build build
cmake --install build
```
The binaries (main program and unit tests) go to ./bin/ inside your build directory.

For a debug build, configure with the ```-DCMAKE_BUILD_TYPE=Debug``` flag.
 

