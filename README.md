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
mkdir build
cd build
cmake ..
make
make install	
```

## Debugging

To build debug version and unit tests use the ```-DCMAKE_BUILD_TYPE=Debug``` flag when calling cmake, ie:

```
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

Google Test framework is used for unit testing. cmake will attempt to automatically fetch the required libraries. 
The binaries (main program and unit tests) go to ./bin/ inside your build directory.
