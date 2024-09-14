# SingletEWPT

C++ program for studying the electroweak phase transition in Standard Model + real singlet theory. Computes the two-loop effective potential in the high-temperature approximation and varies the temperature until phase transition is found. See Arxiv 2103.07467 and 2408.15912 for details on the relevant physics.

The program supports scanning over free parameters of the singlet model. Set parameter ranges in the ```parameters``` file, or use one of the utility scripts in the ```examples/``` folder.

### Requirements

Requires C++17 compliant compiler and CMake 3.15 or newer. We depend on GNU Scientific Library (GSL, https://www.gnu.org/software/gsl/) and ```dlib``` (http://dlib.net). You will need to manually install the GSL libraries. For ```dlib``` you can either manually install the libraries, or let CMake automatically fetch it for you during build stage (see below).

### Compiling and running

```
cmake -B build
cmake --build build
cmake --install build
```
The binaries get installed to ```bin/``` in the project root directory. You need to have ```parameters``` file in the working directory when running the program; this sets parameter ranges and other configuration settings for the program. An example ```parameters``` file is included in repository root. Output of phase transition data goes to ```transitions.dat``` in the working directory and results of zero-temperature checks go to ```data_T0.dat```. Two labels files are also produced that contain description of how to interpret columns of the resulting data. 

For a debug build, configure with the ```-DCMAKE_BUILD_TYPE=Debug``` flag.

Compilation tested with: GCC 9.4.0, GSL 2.5, ```dlib``` 19.24.2.
