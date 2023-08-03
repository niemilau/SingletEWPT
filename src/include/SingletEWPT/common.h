#ifndef COMMON_H
#define COMMON_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <complex>


#define PI 3.14159265358979323846
#define EULERGAMMA 0.57721566490153286060
#define GLAISHER 1.28242712910062263687

// Debug printing macro
#ifdef DEBUG_MODE
#define DEBUG(x) std::cout << x << std::endl;
#else
#define DEBUG(x)
#endif



// Map for storing string, double pairs. Used for practically all action parameter storage and access
using ParameterMap = std::map<std::string, double>;

// Enum for specifying loop order used in various calculations
enum class ELoopOrder {
	tree = 0, loop1 = 1, loop2 = 2
};

// 'Loop order' options for DR. Here it's more complicated than 'tree', '1loop', '2loop' etc
enum class ELoopOrderDR {
    // LO = 1-loop thermal masses only
    LO,
    // NLO = 2-loop thermal masses and 1-loop corrections to couplings
    NLO,
    // NLONo2Loop = same as NLO but with thermal masses at 1-loop only
    NLONo2Loop
};


// Print error message to stderr and exit with an error code
inline void Die(std::string errorMessage, int errorCode) {
	std::cerr << errorMessage;
	exit(errorCode);
}

/* Print a description of all supported options */
inline void PrintUsage(FILE *fp, const char *path) {
    /* take only the last portion of the path */
    const char *basename = strrchr(path, '/');
    basename = basename ? basename + 1 : path;

    fprintf (fp, "usage: %s\n", basename);
}

// check if a string can be interpreted as a number (integer or floating point). Works also with scientific notation
inline bool IsNumber(const std::string string) {
	char* p;
	double converted = strtod(string.c_str(), &p);
    (void)converted; // suppress unused variable warning
	if (*p) {
		// Conversion failed because the input wasn't a number
		return false;
	} else return true;
}

// Safely access an element from map 
template <typename T>
T GetFromMap(const std::map<std::string, T> &map, const std::string &findMe) {

    T value;
    try {
        value = map.at(findMe);
    }
    catch (const std::out_of_range&) {
        Die("Key \"" + findMe + "\" not found in map\n", 432);
    }
    return value;
}


// Print a simple str, double map to stdout
inline void PrintMap(const ParameterMap &map) {

    for(const auto& elem : map) {
		std::cout << elem.first << " " << elem.second << "\n";
	}
}


inline bool FileExists(const std::string& fname) {
    std::ifstream f(fname.c_str());
    return f.good();
}

#endif
