#ifndef TEXT_READER_H
#define TEXT_READER_H

#include "backends/matrix_backends/matrix_backend.h"
#include "parser.h"

namespace SoftComp {

template<typename T=real>
matrix<T> loadtxt(const string &filename, character delim=' ') {
    //! Read a binary file to an array

    std::ifstream file;
    file.open(filename.c_str());

    if(!file) {warn("UNABLE TO READ FILE");}

    integer rows=0, cols;
    std::vector<T> ss;
    for (std::string line; std::getline(file, line); )
    {
        auto words = split(trim(line),delim);
        for (auto &word: words) {
            ss.push_back(str2num<T>(word));
        }
        rows++;
    }
    cols = ss.size()/rows;

    matrix<T> out(rows,cols);

#if defined(HAS_ARMA)
    std::copy(ss.begin(),ss.end(),out.memptr());
#else
    std::copy(ss.begin(),ss.end(),out.data());
#endif
    file.close();
    return out;
}



template<typename T>
boolean savetxt(const matrix<T> &arr, const string &filename, character delim=' ') {
    //! Save an array to a binary file

    std::ofstream file;
    file.open(filename.c_str());

    if(!file) {warn("NO SUCH FILE OR DIRECTORY"); return false;}

    for (auto i=0; i<arr.rows(); ++i) {
        for (auto j=0; j<arr.cols(); ++j) {
            file << arr(i,j) << delim;
        }
        file << "\n";
    }

    return true;
}

}

#endif // TEXT_READER_H
