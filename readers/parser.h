#ifndef PARSER_H
#define PARSER_H

#include "commons/commons.h"
#include <fstream>
#include <sstream>
#include <cctype>
#include <locale>

namespace SoftComp {

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}



void split(const std::string &s, character delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) {
            elems.push_back(item);
        }
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


template<typename T>
T str2num(const std::string &str);
template<>
double str2num<double>(const std::string &str) {
    return std::stod(str);
}
template<>
float str2num<float>(const std::string &str) {
    return std::stof(str);
}
template<>
int str2num<int>(const std::string &str) {
    return std::stoi(str);
}
template<>
long str2num<long>(const std::string &str) {
    return std::stol(str);
}
template<>
long long str2num<long long>(const std::string &str) {
    return std::stoll(str);
}

}

#endif // PARSER_H
