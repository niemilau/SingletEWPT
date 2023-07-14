# SingletEWPT

Program that does stuff

# Dependencies

- dlib C++ library
- Boost odeint

# Installation

```
mkdir build
cd build
cmake ..
make
make install	
```

The 

# Debugging

To build debug version and unit tests use the ```-DCMAKE_BUILD_TYPE=Debug``` flag when calling cmake, ie:

```
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

Google Test framework is used for unit testing. cmake will attempt to automatically fetch the required libraries. 
The binaries (main program and unit tests) go to ./bin/ inside your build directory.
