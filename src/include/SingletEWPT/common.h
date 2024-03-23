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
#include <filesystem>

#include <cctype>    // for std::isspace


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
    
    std::filesystem::path path{ fname };
    return std::filesystem::exists(path);
}

// Remove tabs and whitespaces from start and end of string
inline std::string TrimString(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    size_t last = str.find_last_not_of(" \t\r\n");
    if (first == std::string::npos || last == std::string::npos) {
        // String contains only whitespace
        return "";
    }
    return str.substr(first, last - first + 1);
}


inline std::vector<std::string> SplitString(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}



/* Reads numbers from file, line by line, and puts them in a vector. 
Will not work for files with more than 1 column */
inline std::vector<double> ReadArrayFromFile(const std::string& fileName) {
    std::vector<double> res;
    std::ifstream file(fileName);
    
    res.reserve(100);
    
    if (!file.is_open()) {
        std::cout << "!! Error opening file " << fileName << std::endl;
        return res;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double value;
        
        if (!(iss >> value)) {
            std::cerr << "!! Error parsing line " << line << std::endl;
            continue;
        }
        
        res.push_back(value);
    }
    
    file.close();
    return res;
}


#endif
