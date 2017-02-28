#ifndef BACKEND_CONSTITUTIVETENSORS_H
#define BACKEND_CONSTITUTIVETENSORS_H

#include "commons/commons.h"
#include "backends/matrix_backends/backend_tensor.h"

namespace SoftComp {


tensor<real,3> piolaH(const tensor<real,3> &sigma_H, const tensor<real,3> &F) {
    //--------------------------------------------------------------
    // Concatenate piolaH
    //--------------------------------------------------------------
    tensor<real,3>  out;
    auto ndim    =  F.dimension(1);
    if (ndim==3){
        out =  cross(sigma_H,F)*0.5;
    }
    else {
        out = tensor<real,3>(size(sigma_H,0),2,2);
        real *__restrict__ out_data = out.data();
        const real *__restrict__ sigma_H_data = sigma_H.data();
        for (auto i=0; i<size(sigma_H,0); ++i) {
            out_data[i*4]   =  sigma_H_data[i*4+3];
            out_data[i*4+1] = -sigma_H_data[i*4+2];
            out_data[i*4+2] = -sigma_H_data[i*4+1];;
            out_data[i*4+3] =  sigma_H_data[i*4];
        }
    }
    return out;
}

tensor<real,3> piolaJ(const tensor<real,1> &sigma_J, const tensor<real,3> &H) {
    //--------------------------------------------------------------
    // This function computes the J term in first Piola:
    // (SigmaJ*H)
    //--------------------------------------------------------------
    auto ngauss   =   H.dimension(0);
    auto ndim     =   H.dimension(1);
    tensor<real,3> out(ngauss,ndim,ndim);
    for (integer igauss=0; igauss<ngauss; ++igauss) {
        out.chip(igauss,0)  = sigma_J(igauss)*static_cast<tensor<real,2>>(H.chip(igauss,0));
    }
    return out;
}

//template<typename T>
//SC_INLINE void __Hessian_FH_contribution__(const T *__restrict__ F_data, const T *__restrict__ H_data, T *__restrict__ out_data, const int &ndim) {
////SC_INLINE void __Hessian_FH_contribution__(const tensor<real,3> *__restrict__ F_data, const tensor<real,5> *__restrict__ H_data, tensor<real,5> *__restrict__ out_data, const int &ndim) {
//    //--------------------------------------------------------------
//    // Contribution FH in the Hessian operator
//    //--------------------------------------------------------------
//    if (ndim==3){
//       T F00        =  F_data[0];
//       T F11        =  F_data[4];
//       T F22        =  F_data[8];
//       T F01        =  F_data[1];
//       T F02        =  F_data[2];
//       T F12        =  F_data[5];
//       T F10        =  F_data[3];
//       T F20        =  F_data[6];
//       T F21        =  F_data[7];

//       T W0000  =  H_data[0] ;  T W0001  =  H_data[1] ;  T W0002  =  H_data[2] ;  T W0010  =  H_data[3] ;  T W0011  =  H_data[4] ;
//       T W0012  =  H_data[5] ;  T W0020  =  H_data[6] ;  T W0021  =  H_data[7] ;  T W0022  =  H_data[8] ;  T W0100  =  H_data[9] ;
//       T W0101  =  H_data[10];  T W0102  =  H_data[11];  T W0110  =  H_data[12];  T W0111  =  H_data[13];  T W0112  =  H_data[14];
//       T W0120  =  H_data[15];  T W0121  =  H_data[16];  T W0122  =  H_data[17];  T W0200  =  H_data[18];  T W0201  =  H_data[19];
//       T W0202  =  H_data[20];  T W0210  =  H_data[21];  T W0211  =  H_data[22];  T W0212  =  H_data[23];  T W0220  =  H_data[24];
//       T W0221  =  H_data[25];  T W0222  =  H_data[26];  T W1000  =  H_data[27];  T W1001  =  H_data[28];  T W1002  =  H_data[29];
//       T W1010  =  H_data[30];  T W1011  =  H_data[31];  T W1012  =  H_data[32];  T W1020  =  H_data[33];  T W1021  =  H_data[34];
//       T W1022  =  H_data[35];  T W1100  =  H_data[36];  T W1101  =  H_data[37];  T W1102  =  H_data[38];  T W1110  =  H_data[39];
//       T W1111  =  H_data[40];  T W1112  =  H_data[41];  T W1120  =  H_data[42];  T W1121  =  H_data[43];  T W1122  =  H_data[44];
//       T W1200  =  H_data[45];  T W1201  =  H_data[46];  T W1202  =  H_data[47];  T W1210  =  H_data[48];  T W1211  =  H_data[49];
//       T W1212  =  H_data[50];  T W1220  =  H_data[51];  T W1221  =  H_data[52];  T W1222  =  H_data[53];  T W2000  =  H_data[54];
//       T W2001  =  H_data[55];  T W2002  =  H_data[56];  T W2010  =  H_data[57];  T W2011  =  H_data[58];  T W2012  =  H_data[59];
//       T W2020  =  H_data[60];  T W2021  =  H_data[61];  T W2022  =  H_data[62];  T W2100  =  H_data[63];  T W2101  =  H_data[64];
//       T W2102  =  H_data[65];  T W2110  =  H_data[66];  T W2111  =  H_data[67];  T W2112  =  H_data[68];  T W2120  =  H_data[69];
//       T W2121  =  H_data[70];  T W2122  =  H_data[71];  T W2200  =  H_data[72];  T W2201  =  H_data[73];  T W2202  =  H_data[74];
//       T W2210  =  H_data[75];  T W2211  =  H_data[76];  T W2212  =  H_data[77];  T W2220  =  H_data[78];  T W2221  =  H_data[79];
//       T W2222  =  H_data[80];

//       out_data[0]   =  F11*W0022 - F12*W0021 - F21*W0012 + F22*W0011;
//       out_data[1]   =  F12*W0020 - F10*W0022 + F20*W0012 - F22*W0010;
//       out_data[2]   =  F10*W0021 - F11*W0020 - F20*W0011 + F21*W0010;
//       out_data[3]   =  F02*W0021 - F01*W0022 + F21*W0002 - F22*W0001;
//       out_data[4]   =  F00*W0022 - F02*W0020 - F20*W0002 + F22*W0000;
//       out_data[5]   =  F01*W0020 - F00*W0021 + F20*W0001 - F21*W0000;
//       out_data[6]   =  F01*W0012 - F02*W0011 - F11*W0002 + F12*W0001;
//       out_data[7]   =  F02*W0010 - F00*W0012 + F10*W0002 - F12*W0000;
//       out_data[8]   =  F00*W0011 - F01*W0010 - F10*W0001 + F11*W0000;
//       out_data[9]   =  F11*W0122 - F12*W0121 - F21*W0112 + F22*W0111;
//       out_data[10]  =  F12*W0120 - F10*W0122 + F20*W0112 - F22*W0110;
//       out_data[11]  =  F10*W0121 - F11*W0120 - F20*W0111 + F21*W0110;
//       out_data[12]  =  F02*W0121 - F01*W0122 + F21*W0102 - F22*W0101;
//       out_data[13]  =  F00*W0122 - F02*W0120 - F20*W0102 + F22*W0100;
//       out_data[14]  =  F01*W0120 - F00*W0121 + F20*W0101 - F21*W0100;
//       out_data[15]  =  F01*W0112 - F02*W0111 - F11*W0102 + F12*W0101;
//       out_data[16]  =  F02*W0110 - F00*W0112 + F10*W0102 - F12*W0100;
//       out_data[17]  =  F00*W0111 - F01*W0110 - F10*W0101 + F11*W0100;
//       out_data[18]  =  F11*W0222 - F12*W0221 - F21*W0212 + F22*W0211;
//       out_data[19]  =  F12*W0220 - F10*W0222 + F20*W0212 - F22*W0210;
//       out_data[20]  =  F10*W0221 - F11*W0220 - F20*W0211 + F21*W0210;
//       out_data[21]  =  F02*W0221 - F01*W0222 + F21*W0202 - F22*W0201;
//       out_data[22]  =  F00*W0222 - F02*W0220 - F20*W0202 + F22*W0200;
//       out_data[23]  =  F01*W0220 - F00*W0221 + F20*W0201 - F21*W0200;
//       out_data[24]  =  F01*W0212 - F02*W0211 - F11*W0202 + F12*W0201;
//       out_data[25]  =  F02*W0210 - F00*W0212 + F10*W0202 - F12*W0200;
//       out_data[26]  =  F00*W0211 - F01*W0210 - F10*W0201 + F11*W0200;
//       out_data[27]  =  F11*W1022 - F12*W1021 - F21*W1012 + F22*W1011;
//       out_data[28]  =  F12*W1020 - F10*W1022 + F20*W1012 - F22*W1010;
//       out_data[29]  =  F10*W1021 - F11*W1020 - F20*W1011 + F21*W1010;
//       out_data[30]  =  F02*W1021 - F01*W1022 + F21*W1002 - F22*W1001;
//       out_data[31]  =  F00*W1022 - F02*W1020 - F20*W1002 + F22*W1000;
//       out_data[32]  =  F01*W1020 - F00*W1021 + F20*W1001 - F21*W1000;
//       out_data[33]  =  F01*W1012 - F02*W1011 - F11*W1002 + F12*W1001;
//       out_data[34]  =  F02*W1010 - F00*W1012 + F10*W1002 - F12*W1000;
//       out_data[35]  =  F00*W1011 - F01*W1010 - F10*W1001 + F11*W1000;
//       out_data[36]  =  F11*W1122 - F12*W1121 - F21*W1112 + F22*W1111;
//       out_data[37]  =  F12*W1120 - F10*W1122 + F20*W1112 - F22*W1110;
//       out_data[38]  =  F10*W1121 - F11*W1120 - F20*W1111 + F21*W1110;
//       out_data[39]  =  F02*W1121 - F01*W1122 + F21*W1102 - F22*W1101;
//       out_data[40]  =  F00*W1122 - F02*W1120 - F20*W1102 + F22*W1100;
//       out_data[41]  =  F01*W1120 - F00*W1121 + F20*W1101 - F21*W1100;
//       out_data[42]  =  F01*W1112 - F02*W1111 - F11*W1102 + F12*W1101;
//       out_data[43]  =  F02*W1110 - F00*W1112 + F10*W1102 - F12*W1100;
//       out_data[44]  =  F00*W1111 - F01*W1110 - F10*W1101 + F11*W1100;
//       out_data[45]  =  F11*W1222 - F12*W1221 - F21*W1212 + F22*W1211;
//       out_data[46]  =  F12*W1220 - F10*W1222 + F20*W1212 - F22*W1210;
//       out_data[47]  =  F10*W1221 - F11*W1220 - F20*W1211 + F21*W1210;
//       out_data[48]  =  F02*W1221 - F01*W1222 + F21*W1202 - F22*W1201;
//       out_data[49]  =  F00*W1222 - F02*W1220 - F20*W1202 + F22*W1200;
//       out_data[50]  =  F01*W1220 - F00*W1221 + F20*W1201 - F21*W1200;
//       out_data[51]  =  F01*W1212 - F02*W1211 - F11*W1202 + F12*W1201;
//       out_data[52]  =  F02*W1210 - F00*W1212 + F10*W1202 - F12*W1200;
//       out_data[53]  =  F00*W1211 - F01*W1210 - F10*W1201 + F11*W1200;
//       out_data[54]  =  F11*W2022 - F12*W2021 - F21*W2012 + F22*W2011;
//       out_data[55]  =  F12*W2020 - F10*W2022 + F20*W2012 - F22*W2010;
//       out_data[56]  =  F10*W2021 - F11*W2020 - F20*W2011 + F21*W2010;
//       out_data[57]  =  F02*W2021 - F01*W2022 + F21*W2002 - F22*W2001;
//       out_data[58]  =  F00*W2022 - F02*W2020 - F20*W2002 + F22*W2000;
//       out_data[59]  =  F01*W2020 - F00*W2021 + F20*W2001 - F21*W2000;
//       out_data[60]  =  F01*W2012 - F02*W2011 - F11*W2002 + F12*W2001;
//       out_data[61]  =  F02*W2010 - F00*W2012 + F10*W2002 - F12*W2000;
//       out_data[62]  =  F00*W2011 - F01*W2010 - F10*W2001 + F11*W2000;
//       out_data[63]  =  F11*W2122 - F12*W2121 - F21*W2112 + F22*W2111;
//       out_data[64]  =  F12*W2120 - F10*W2122 + F20*W2112 - F22*W2110;
//       out_data[65]  =  F10*W2121 - F11*W2120 - F20*W2111 + F21*W2110;
//       out_data[66]  =  F02*W2121 - F01*W2122 + F21*W2102 - F22*W2101;
//       out_data[67]  =  F00*W2122 - F02*W2120 - F20*W2102 + F22*W2100;
//       out_data[68]  =  F01*W2120 - F00*W2121 + F20*W2101 - F21*W2100;
//       out_data[69]  =  F01*W2112 - F02*W2111 - F11*W2102 + F12*W2101;
//       out_data[70]  =  F02*W2110 - F00*W2112 + F10*W2102 - F12*W2100;
//       out_data[71]  =  F00*W2111 - F01*W2110 - F10*W2101 + F11*W2100;
//       out_data[72]  =  F11*W2222 - F12*W2221 - F21*W2212 + F22*W2211;
//       out_data[73]  =  F12*W2220 - F10*W2222 + F20*W2212 - F22*W2210;
//       out_data[74]  =  F10*W2221 - F11*W2220 - F20*W2211 + F21*W2210;
//       out_data[75]  =  F02*W2221 - F01*W2222 + F21*W2202 - F22*W2201;
//       out_data[76]  =  F00*W2222 - F02*W2220 - F20*W2202 + F22*W2200;
//       out_data[77]  =  F01*W2220 - F00*W2221 + F20*W2201 - F21*W2200;
//       out_data[78]  =  F01*W2212 - F02*W2211 - F11*W2202 + F12*W2201;
//       out_data[79]  =  F02*W2210 - F00*W2212 + F10*W2202 - F12*W2200;
//       out_data[80]  =  F00*W2211 - F01*W2210 - F10*W2201 + F11*W2200;
//    }
//    else{
//       T W0000  =  H_data[0];   T W0001  =  H_data[1];   T W0010  =  H_data[2];   T W0011  =  H_data[3];   T W0100  =  H_data[4];
//       T W0101  =  H_data[5];   T W0110  =  H_data[6];   T W0111  =  H_data[7];   T W1000  =  H_data[8];   T W1001  =  H_data[9];
//       T W1010  =  H_data[10];  T W1011  =  H_data[11];  T W1100  =  H_data[12];  T W1101  =  H_data[13];  T W1110  =  H_data[14];
//       T W1111  =  H_data[15];

//       out_data[0]   =  W0011;   out_data[1]   =  -W0010;   out_data[2]   =  -W0001;   out_data[3]  =  W0000;
//       out_data[4]   =  W0111;   out_data[5]   =  -W0110;   out_data[6]   =  -W0101;   out_data[7]  =  W0100;
//       out_data[8]   =  W1011;   out_data[9]   =  -W1010;   out_data[10]  =  -W1001;  out_data[11]  =  W1000;
//       out_data[12]  =  W1111;  out_data[13]   =  -W1110;  out_data[14]   =  -W1101;  out_data[15]  =  W1100;
//    }
//}

//tensor<real,3> H_FH(const tensor<real,3> &F, const tensor<real,5> &WFH,  const int &ndim) {
//    //--------------------------------------------------------------
//    // Obtain FH contribution in the Hessian for every Gauss point
//    //--------------------------------------------------------------
//    auto ngauss        =  F.dimension(0);
//    tensor<real,5> out(ngauss,ndim,ndim,ndim,ndim);
//    for (auto i=0; i<ngauss; ++i) {
//        tensor<real,2> F_2d   =  F.chip(i,0);
//        tensor<real,2> WFH_4d =  WFH.chip(i,0);
//        __Hessian_FH_contribution__(F_2d.data(),WFH_4d.data(), out.data()+ndim*ndim*ndim*ndim*i,ndim);
//    }
//    return out;
//}

//template<typename T>
//SC_INLINE void __Hessian_FJ_contribution__(const T *__restrict__ a_data, const T *__restrict__ b_data, T *__restrict__ out_data, const int &ndim) {
//    //--------------------------------------------------------------
//    // Contribution FJ in the Hessian operator
//    //--------------------------------------------------------------
//    if (ndim==3){
//       T H00         =  a_data[0];
//       T H11         =  a_data[4];
//       T H22         =  a_data[8];
//       T H01         =  a_data[1];
//       T H02         =  a_data[2];
//       T H12         =  a_data[5];
//       T H10         =  a_data[3];
//       T H20         =  a_data[6];
//       T H21         =  a_data[7];

//       T WFJ00       =  b_data[0];
//       T WFJ11       =  b_data[4];
//       T WFJ22       =  b_data[8];
//       T WFJ01       =  b_data[1];
//       T WFJ02       =  b_data[2];
//       T WFJ12       =  b_data[5];
//       T WFJ10       =  b_data[3];
//       T WFJ20       =  b_data[6];
//       T WFJ21       =  b_data[7];

//       out_data[0]   =  H00*WFJ00;
//       out_data[1]   =  H01*WFJ00;
//       out_data[2]   =  H02*WFJ00;
//       out_data[3]   =  H10*WFJ00;
//       out_data[4]   =  H11*WFJ00;
//       out_data[5]   =  H12*WFJ00;
//       out_data[6]   =  H20*WFJ00;
//       out_data[7]   =  H21*WFJ00;
//       out_data[8]   =  H22*WFJ00;
//       out_data[9]   =  H00*WFJ01;
//       out_data[10]  =  H01*WFJ01;
//       out_data[11]  =  H02*WFJ01;
//       out_data[12]  =  H10*WFJ01;
//       out_data[13]  =  H11*WFJ01;
//       out_data[14]  =  H12*WFJ01;
//       out_data[15]  =  H20*WFJ01;
//       out_data[16]  =  H21*WFJ01;
//       out_data[17]  =  H22*WFJ01;
//       out_data[18]  =  H00*WFJ02;
//       out_data[19]  =  H01*WFJ02;
//       out_data[20]  =  H02*WFJ02;
//       out_data[21]  =  H10*WFJ02;
//       out_data[22]  =  H11*WFJ02;
//       out_data[23]  =  H12*WFJ02;
//       out_data[24]  =  H20*WFJ02;
//       out_data[25]  =  H21*WFJ02;
//       out_data[26]  =  H22*WFJ02;
//       out_data[27]  =  H00*WFJ10;
//       out_data[28]  =  H01*WFJ10;
//       out_data[29]  =  H02*WFJ10;
//       out_data[30]  =  H10*WFJ10;
//       out_data[31]  =  H11*WFJ10;
//       out_data[32]  =  H12*WFJ10;
//       out_data[33]  =  H20*WFJ10;
//       out_data[34]  =  H21*WFJ10;
//       out_data[35]  =  H22*WFJ10;
//       out_data[36]  =  H00*WFJ11;
//       out_data[37]  =  H01*WFJ11;
//       out_data[38]  =  H02*WFJ11;
//       out_data[39]  =  H10*WFJ11;
//       out_data[40]  =  H11*WFJ11;
//       out_data[41]  =  H12*WFJ11;
//       out_data[42]  =  H20*WFJ11;
//       out_data[43]  =  H21*WFJ11;
//       out_data[44]  =  H22*WFJ11;
//       out_data[45]  =  H00*WFJ12;
//       out_data[46]  =  H01*WFJ12;
//       out_data[47]  =  H02*WFJ12;
//       out_data[48]  =  H10*WFJ12;
//       out_data[49]  =  H11*WFJ12;
//       out_data[50]  =  H12*WFJ12;
//       out_data[51]  =  H20*WFJ12;
//       out_data[52]  =  H21*WFJ12;
//       out_data[53]  =  H22*WFJ12;
//       out_data[54]  =  H00*WFJ20;
//       out_data[55]  =  H01*WFJ20;
//       out_data[56]  =  H02*WFJ20;
//       out_data[57]  =  H10*WFJ20;
//       out_data[58]  =  H11*WFJ20;
//       out_data[59]  =  H12*WFJ20;
//       out_data[60]  =  H20*WFJ20;
//       out_data[61]  =  H21*WFJ20;
//       out_data[62]  =  H22*WFJ20;
//       out_data[63]  =  H00*WFJ21;
//       out_data[64]  =  H01*WFJ21;
//       out_data[65]  =  H02*WFJ21;
//       out_data[66]  =  H10*WFJ21;
//       out_data[67]  =  H11*WFJ21;
//       out_data[68]  =  H12*WFJ21;
//       out_data[69]  =  H20*WFJ21;
//       out_data[70]  =  H21*WFJ21;
//       out_data[71]  =  H22*WFJ21;
//       out_data[72]  =  H00*WFJ22;
//       out_data[73]  =  H01*WFJ22;
//       out_data[74]  =  H02*WFJ22;
//       out_data[75]  =  H10*WFJ22;
//       out_data[76]  =  H11*WFJ22;
//       out_data[77]  =  H12*WFJ22;
//       out_data[78]  =  H20*WFJ22;
//       out_data[79]  =  H21*WFJ22;
//       out_data[80]  =  H22*WFJ22;
//    }
//    else{
//        T H00         =  a_data[0];
//        T H11         =  a_data[3];
//        T H01         =  a_data[1];
//        T H10         =  a_data[2];

//        T WFJ00       =  b_data[0];
//        T WFJ11       =  b_data[3];
//        T WFJ01       =  b_data[1];
//        T WFJ10       =  b_data[2];

//        out_data[0]   =  H00*WFJ00;
//        out_data[1]   =  H01*WFJ00;
//        out_data[2]   =  H10*WFJ00;
//        out_data[3]   =  H11*WFJ00;
//        out_data[4]   =  H00*WFJ01;
//        out_data[5]   =  H01*WFJ01;
//        out_data[6]   =  H10*WFJ01;
//        out_data[7]   =  H11*WFJ01;
//        out_data[8]   =  H00*WFJ10;
//        out_data[9]   =  H01*WFJ10;
//        out_data[10]  =  H10*WFJ10;
//        out_data[11]  =  H11*WFJ10;
//        out_data[12]  =  H00*WFJ11;
//        out_data[13]  =  H01*WFJ11;
//        out_data[14]  =  H10*WFJ11;
//        out_data[15]  =  H11*WFJ11;
//    }
//}

//tensor<real,3> H_FJ(const tensor<real,3> &H, const tensor<real,3> &WFJ,  const int &ndim) {
//    //--------------------------------------------------------------
//    // Obtain FJ contribution in the Hessian for every Gauss point
//    //--------------------------------------------------------------
//    auto ngauss        =  H.dimension(2);
//    tensor<real,5> out(ndim,ndim,ndim,ndim,ngauss);
//    for (auto i=0; i<ngauss; ++i) {
//        tensor<real,2> H_2d   =  H.chip(i,2);
//        tensor<real,2> WFJ_2d =  WFJ.chip(i,2);
//        __Hessian_FJ_contribution__(H_2d.data(),WFJ.data(), out.data()+ndim*ndim*ndim*ndim*i,ndim);
//    }
//    return out;
//}

//template<typename T>
//SC_INLINE void __Hessian_HF_contribution__(const T *__restrict__ F_data, const T *__restrict__ H_data, T *__restrict__ out_data, const int &ndim) {
//    //--------------------------------------------------------------
//    // Contribution HF in the Hessian operator
//    //--------------------------------------------------------------
//    if (ndim==3){
//       T F00        =  F_data[0];
//       T F11        =  F_data[4];
//       T F22        =  F_data[8];
//       T F01        =  F_data[1];
//       T F02        =  F_data[2];
//       T F12        =  F_data[5];
//       T F10        =  F_data[3];
//       T F20        =  F_data[6];
//       T F21        =  F_data[7];

//       T W0000  =  H_data[0] ;  T W0001  =  H_data[1] ;  T W0002  =  H_data[2] ;  T W0010  =  H_data[3] ;  T W0011  =  H_data[4] ;
//       T W0012  =  H_data[5] ;  T W0020  =  H_data[6] ;  T W0021  =  H_data[7] ;  T W0022  =  H_data[8] ;  T W0100  =  H_data[9] ;
//       T W0101  =  H_data[10];  T W0102  =  H_data[11];  T W0110  =  H_data[12];  T W0111  =  H_data[13];  T W0112  =  H_data[14];
//       T W0120  =  H_data[15];  T W0121  =  H_data[16];  T W0122  =  H_data[17];  T W0200  =  H_data[18];  T W0201  =  H_data[19];
//       T W0202  =  H_data[20];  T W0210  =  H_data[21];  T W0211  =  H_data[22];  T W0212  =  H_data[23];  T W0220  =  H_data[24];
//       T W0221  =  H_data[25];  T W0222  =  H_data[26];  T W1000  =  H_data[27];  T W1001  =  H_data[28];  T W1002  =  H_data[29];
//       T W1010  =  H_data[30];  T W1011  =  H_data[31];  T W1012  =  H_data[32];  T W1020  =  H_data[33];  T W1021  =  H_data[34];
//       T W1022  =  H_data[35];  T W1100  =  H_data[36];  T W1101  =  H_data[37];  T W1102  =  H_data[38];  T W1110  =  H_data[39];
//       T W1111  =  H_data[40];  T W1112  =  H_data[41];  T W1120  =  H_data[42];  T W1121  =  H_data[43];  T W1122  =  H_data[44];
//       T W1200  =  H_data[45];  T W1201  =  H_data[46];  T W1202  =  H_data[47];  T W1210  =  H_data[48];  T W1211  =  H_data[49];
//       T W1212  =  H_data[50];  T W1220  =  H_data[51];  T W1221  =  H_data[52];  T W1222  =  H_data[53];  T W2000  =  H_data[54];
//       T W2001  =  H_data[55];  T W2002  =  H_data[56];  T W2010  =  H_data[57];  T W2011  =  H_data[58];  T W2012  =  H_data[59];
//       T W2020  =  H_data[60];  T W2021  =  H_data[61];  T W2022  =  H_data[62];  T W2100  =  H_data[63];  T W2101  =  H_data[64];
//       T W2102  =  H_data[65];  T W2110  =  H_data[66];  T W2111  =  H_data[67];  T W2112  =  H_data[68];  T W2120  =  H_data[69];
//       T W2121  =  H_data[70];  T W2122  =  H_data[71];  T W2200  =  H_data[72];  T W2201  =  H_data[73];  T W2202  =  H_data[74];
//       T W2210  =  H_data[75];  T W2211  =  H_data[76];  T W2212  =  H_data[77];  T W2220  =  H_data[78];  T W2221  =  H_data[79];
//       T W2222  =  H_data[80];

//       out_data[0]   =  F22*W1100 - F21*W1200 - F12*W2100 + F11*W2200;
//       out_data[1]   =  F22*W1101 - F21*W1201 - F12*W2101 + F11*W2201;
//       out_data[2]   =  F22*W1102 - F21*W1202 - F12*W2102 + F11*W2202;
//       out_data[3]   =  F22*W1110 - F21*W1210 - F12*W2110 + F11*W2210;
//       out_data[4]   =  F22*W1111 - F21*W1211 - F12*W2111 + F11*W2211;
//       out_data[5]   =  F22*W1112 - F21*W1212 - F12*W2112 + F11*W2212;
//       out_data[6]   =  F22*W1120 - F21*W1220 - F12*W2120 + F11*W2220;
//       out_data[7]   =  F22*W1121 - F21*W1221 - F12*W2121 + F11*W2221;
//       out_data[8]   =  F22*W1122 - F21*W1222 - F12*W2122 + F11*W2222;
//       out_data[9]   =  F20*W1200 - F22*W1000 + F12*W2000 - F10*W2200;
//       out_data[10]  =  F20*W1201 - F22*W1001 + F12*W2001 - F10*W2201;
//       out_data[11]  =  F20*W1202 - F22*W1002 + F12*W2002 - F10*W2202;
//       out_data[12]  =  F20*W1210 - F22*W1010 + F12*W2010 - F10*W2210;
//       out_data[13]  =  F20*W1211 - F22*W1011 + F12*W2011 - F10*W2211;
//       out_data[14]  =  F20*W1212 - F22*W1012 + F12*W2012 - F10*W2212;
//       out_data[15]  =  F20*W1220 - F22*W1020 + F12*W2020 - F10*W2220;
//       out_data[16]  =  F20*W1221 - F22*W1021 + F12*W2021 - F10*W2221;
//       out_data[17]  =  F20*W1222 - F22*W1022 + F12*W2022 - F10*W2222;
//       out_data[18]  =  F21*W1000 - F20*W1100 - F11*W2000 + F10*W2100;
//       out_data[19]  =  F21*W1001 - F20*W1101 - F11*W2001 + F10*W2101;
//       out_data[20]  =  F21*W1002 - F20*W1102 - F11*W2002 + F10*W2102;
//       out_data[21]  =  F21*W1010 - F20*W1110 - F11*W2010 + F10*W2110;
//       out_data[22]  =  F21*W1011 - F20*W1111 - F11*W2011 + F10*W2111;
//       out_data[23]  =  F21*W1012 - F20*W1112 - F11*W2012 + F10*W2112;
//       out_data[24]  =  F21*W1020 - F20*W1120 - F11*W2020 + F10*W2120;
//       out_data[25]  =  F21*W1021 - F20*W1121 - F11*W2021 + F10*W2121;
//       out_data[26]  =  F21*W1022 - F20*W1122 - F11*W2022 + F10*W2122;
//       out_data[27]  =  F21*W0200 - F22*W0100 + F02*W2100 - F01*W2200;
//       out_data[28]  =  F21*W0201 - F22*W0101 + F02*W2101 - F01*W2201;
//       out_data[29]  =  F21*W0202 - F22*W0102 + F02*W2102 - F01*W2202;
//       out_data[30]  =  F21*W0210 - F22*W0110 + F02*W2110 - F01*W2210;
//       out_data[31]  =  F21*W0211 - F22*W0111 + F02*W2111 - F01*W2211;
//       out_data[32]  =  F21*W0212 - F22*W0112 + F02*W2112 - F01*W2212;
//       out_data[33]  =  F21*W0220 - F22*W0120 + F02*W2120 - F01*W2220;
//       out_data[34]  =  F21*W0221 - F22*W0121 + F02*W2121 - F01*W2221;
//       out_data[35]  =  F21*W0222 - F22*W0122 + F02*W2122 - F01*W2222;
//       out_data[36]  =  F22*W0000 - F20*W0200 - F02*W2000 + F00*W2200;
//       out_data[37]  =  F22*W0001 - F20*W0201 - F02*W2001 + F00*W2201;
//       out_data[38]  =  F22*W0002 - F20*W0202 - F02*W2002 + F00*W2202;
//       out_data[39]  =  F22*W0010 - F20*W0210 - F02*W2010 + F00*W2210;
//       out_data[40]  =  F22*W0011 - F20*W0211 - F02*W2011 + F00*W2211;
//       out_data[41]  =  F22*W0012 - F20*W0212 - F02*W2012 + F00*W2212;
//       out_data[42]  =  F22*W0020 - F20*W0220 - F02*W2020 + F00*W2220;
//       out_data[43]  =  F22*W0021 - F20*W0221 - F02*W2021 + F00*W2221;
//       out_data[44]  =  F22*W0022 - F20*W0222 - F02*W2022 + F00*W2222;
//       out_data[45]  =  F20*W0100 - F21*W0000 + F01*W2000 - F00*W2100;
//       out_data[46]  =  F20*W0101 - F21*W0001 + F01*W2001 - F00*W2101;
//       out_data[47]  =  F20*W0102 - F21*W0002 + F01*W2002 - F00*W2102;
//       out_data[48]  =  F20*W0110 - F21*W0010 + F01*W2010 - F00*W2110;
//       out_data[49]  =  F20*W0111 - F21*W0011 + F01*W2011 - F00*W2111;
//       out_data[50]  =  F20*W0112 - F21*W0012 + F01*W2012 - F00*W2112;
//       out_data[51]  =  F20*W0120 - F21*W0020 + F01*W2020 - F00*W2120;
//       out_data[52]  =  F20*W0121 - F21*W0021 + F01*W2021 - F00*W2121;
//       out_data[53]  =  F20*W0122 - F21*W0022 + F01*W2022 - F00*W2122;
//       out_data[54]  =  F12*W0100 - F11*W0200 - F02*W1100 + F01*W1200;
//       out_data[55]  =  F12*W0101 - F11*W0201 - F02*W1101 + F01*W1201;
//       out_data[56]  =  F12*W0102 - F11*W0202 - F02*W1102 + F01*W1202;
//       out_data[57]  =  F12*W0110 - F11*W0210 - F02*W1110 + F01*W1210;
//       out_data[58]  =  F12*W0111 - F11*W0211 - F02*W1111 + F01*W1211;
//       out_data[59]  =  F12*W0112 - F11*W0212 - F02*W1112 + F01*W1212;
//       out_data[60]  =  F12*W0120 - F11*W0220 - F02*W1120 + F01*W1220;
//       out_data[61]  =  F12*W0121 - F11*W0221 - F02*W1121 + F01*W1221;
//       out_data[62]  =  F12*W0122 - F11*W0222 - F02*W1122 + F01*W1222;
//       out_data[63]  =  F10*W0200 - F12*W0000 + F02*W1000 - F00*W1200;
//       out_data[64]  =  F10*W0201 - F12*W0001 + F02*W1001 - F00*W1201;
//       out_data[65]  =  F10*W0202 - F12*W0002 + F02*W1002 - F00*W1202;
//       out_data[66]  =  F10*W0210 - F12*W0010 + F02*W1010 - F00*W1210;
//       out_data[67]  =  F10*W0211 - F12*W0011 + F02*W1011 - F00*W1211;
//       out_data[68]  =  F10*W0212 - F12*W0012 + F02*W1012 - F00*W1212;
//       out_data[69]  =  F10*W0220 - F12*W0020 + F02*W1020 - F00*W1220;
//       out_data[70]  =  F10*W0221 - F12*W0021 + F02*W1021 - F00*W1221;
//       out_data[71]  =  F10*W0222 - F12*W0022 + F02*W1022 - F00*W1222;
//       out_data[72]  =  F11*W0000 - F10*W0100 - F01*W1000 + F00*W1100;
//       out_data[73]  =  F11*W0001 - F10*W0101 - F01*W1001 + F00*W1101;
//       out_data[74]  =  F11*W0002 - F10*W0102 - F01*W1002 + F00*W1102;
//       out_data[75]  =  F11*W0010 - F10*W0110 - F01*W1010 + F00*W1110;
//       out_data[76]  =  F11*W0011 - F10*W0111 - F01*W1011 + F00*W1111;
//       out_data[77]  =  F11*W0012 - F10*W0112 - F01*W1012 + F00*W1112;
//       out_data[78]  =  F11*W0020 - F10*W0120 - F01*W1020 + F00*W1120;
//       out_data[79]  =  F11*W0021 - F10*W0121 - F01*W1021 + F00*W1121;
//       out_data[80]  =  F11*W0022 - F10*W0122 - F01*W1022 + F00*W1122;
//    }
//    else{
//       T W0000  =  H_data[0];   T W0001  =  H_data[1];   T W0010  =  H_data[2];   T W0011  =  H_data[3];   T W0100  =  H_data[4];
//       T W0101  =  H_data[5];   T W0110  =  H_data[6];   T W0111  =  H_data[7];   T W1000  =  H_data[8];   T W1001  =  H_data[9];
//       T W1010  =  H_data[10];  T W1011  =  H_data[11];  T W1100  =  H_data[12];  T W1101  =  H_data[13];  T W1110  =  H_data[14];
//       T W1111  =  H_data[15];

//       out_data[0]   =   W1100;   out_data[1]   =   W1101;   out_data[2]   =   W1110;   out_data[3]   =   W1111;
//       out_data[4]   =  -W1000;   out_data[5]   =  -W1001;   out_data[6]   =  -W1010;   out_data[7]   =  -W1011;
//       out_data[8]   =  -W0100;   out_data[9]   =  -W0101;   out_data[10]  =  -W0110;   out_data[11]  =  -W0111;
//       out_data[12]  =   W0000;   out_data[13]  =   W0001;   out_data[14]  =   W0010;   out_data[15]  =   W0011;
//    }
//}

//tensor<real,3> H_HF(const tensor<real,3> &F, const tensor<real,5> &WHF,  const int &ndim) {
//    //--------------------------------------------------------------
//    // Obtain HF contribution in the Hessian for every Gauss point
//    //--------------------------------------------------------------
//    auto ngauss        =  F.dimension(2);
//    tensor<real,5> out(ndim,ndim,ndim,ndim,ngauss);
//    for (auto i=0; i<ngauss; ++i) {
//        tensor<real,2> F_2d   =  F.chip(i,2);
//        tensor<real,2> WHF_4d =  WHF.chip(i,4);
//        __Hessian_HF_contribution__(F_2d.data(),WHF_4d.data(), out.data()+ndim*ndim*ndim*ndim*i,ndim);
//    }
//    return out;
//}

template<typename T>
SC_INLINE void __Hessian_HH_contribution__(const T *__restrict__ F_data, const T *__restrict__ H_data, T *__restrict__ out_data, const int &ndim) {
    //--------------------------------------------------------------
    // Contribution HF in the Hessian operator
    //--------------------------------------------------------------
    if (ndim==3){
       T F00        =  F_data[0];
       T F11        =  F_data[4];
       T F22        =  F_data[8];
       T F01        =  F_data[1];
       T F02        =  F_data[2];
       T F12        =  F_data[5];
       T F10        =  F_data[3];
       T F20        =  F_data[6];
       T F21        =  F_data[7];

       T W0000  =  H_data[0] ;  T W0001  =  H_data[1] ;  T W0002  =  H_data[2] ;  T W0010  =  H_data[3] ;  T W0011  =  H_data[4] ;
       T W0012  =  H_data[5] ;  T W0020  =  H_data[6] ;  T W0021  =  H_data[7] ;  T W0022  =  H_data[8] ;  T W0100  =  H_data[9] ;
       T W0101  =  H_data[10];  T W0102  =  H_data[11];  T W0110  =  H_data[12];  T W0111  =  H_data[13];  T W0112  =  H_data[14];
       T W0120  =  H_data[15];  T W0121  =  H_data[16];  T W0122  =  H_data[17];  T W0200  =  H_data[18];  T W0201  =  H_data[19];
       T W0202  =  H_data[20];  T W0210  =  H_data[21];  T W0211  =  H_data[22];  T W0212  =  H_data[23];  T W0220  =  H_data[24];
       T W0221  =  H_data[25];  T W0222  =  H_data[26];  T W1000  =  H_data[27];  T W1001  =  H_data[28];  T W1002  =  H_data[29];
       T W1010  =  H_data[30];  T W1011  =  H_data[31];  T W1012  =  H_data[32];  T W1020  =  H_data[33];  T W1021  =  H_data[34];
       T W1022  =  H_data[35];  T W1100  =  H_data[36];  T W1101  =  H_data[37];  T W1102  =  H_data[38];  T W1110  =  H_data[39];
       T W1111  =  H_data[40];  T W1112  =  H_data[41];  T W1120  =  H_data[42];  T W1121  =  H_data[43];  T W1122  =  H_data[44];
       T W1200  =  H_data[45];  T W1201  =  H_data[46];  T W1202  =  H_data[47];  T W1210  =  H_data[48];  T W1211  =  H_data[49];
       T W1212  =  H_data[50];  T W1220  =  H_data[51];  T W1221  =  H_data[52];  T W1222  =  H_data[53];  T W2000  =  H_data[54];
       T W2001  =  H_data[55];  T W2002  =  H_data[56];  T W2010  =  H_data[57];  T W2011  =  H_data[58];  T W2012  =  H_data[59];
       T W2020  =  H_data[60];  T W2021  =  H_data[61];  T W2022  =  H_data[62];  T W2100  =  H_data[63];  T W2101  =  H_data[64];
       T W2102  =  H_data[65];  T W2110  =  H_data[66];  T W2111  =  H_data[67];  T W2112  =  H_data[68];  T W2120  =  H_data[69];
       T W2121  =  H_data[70];  T W2122  =  H_data[71];  T W2200  =  H_data[72];  T W2201  =  H_data[73];  T W2202  =  H_data[74];
       T W2210  =  H_data[75];  T W2211  =  H_data[76];  T W2212  =  H_data[77];  T W2220  =  H_data[78];  T W2221  =  H_data[79];
       T W2222  =  H_data[80];

       //out_data[0]   =  F22***W1111 + F21^2*W1212 +   F12^2*W2121 + F11^2*W2222 + F11*F22*W1122 - F12*F22*W1121 - F21*F22*W1112 - F11*F21*W1222 + F12*F21*W1221 - F21*F22*W1211 - F11*F12*W2122 + F12*F21*W2112 - F12*F22*W2111 - F11*F12*W2221 - F11*F21*W2212 + F11*F22*W2211;
       out_data[0]   =  F22*F22*W1111 + F21*F21*W1212 + F12*F12*W2121 + F11*F11*W2222 + F11*F22*W1122 - F12*F22*W1121 - F21*F22*W1112 - F11*F21*W1222 + F12*F21*W1221 - F21*F22*W1211 - F11*F12*W2122 + F12*F21*W2112 - F12*F22*W2111 - F11*F12*W2221 - F11*F21*W2212 + F11*F22*W2211;
       out_data[1]   =  F12*F22*W1120 - F12*F12*W2120 - F10*F22*W1122 - F22*F22*W1110 + F20*F22*W1112 + F10*F21*W1222 - F12*F21*W1220 - F20*F21*W1212 + F21*F22*W1210 + F10*F12*W2122 - F12*F20*W2112 + F12*F22*W2110 - F10*F11*W2222 + F11*F12*W2220 + F11*F20*W2212 - F11*F22*W2210;
       out_data[2]   =  F10*F22*W1121 - F11*F11*W2220 - F21*F21*W1210 - F11*F22*W1120 - F20*F22*W1111 + F21*F22*W1110 - F10*F21*W1221 + F11*F21*W1220 + F20*F21*W1211 - F10*F12*W2121 + F11*F12*W2120 + F12*F20*W2111 - F12*F21*W2110 + F10*F11*W2221 - F11*F20*W2211 + F11*F21*W2210;
       out_data[3]   =  F02*F22*W1121 - F21*F21*W1202 - F01*F22*W1122 - F22*F22*W1101 + F21*F22*W1102 + F01*F21*W1222 - F02*F21*W1221 + F21*F22*W1201 + F01*F12*W2122 - F02*F12*W2121 - F12*F21*W2102 + F12*F22*W2101 - F01*F11*W2222 + F02*F11*W2221 + F11*F21*W2202 - F11*F22*W2201;
       out_data[4]   =  F22*F22*W1100 + F00*F22*W1122 - F02*F22*W1120 - F20*F22*W1102 - F00*F21*W1222 + F02*F21*W1220 + F20*F21*W1202 - F21*F22*W1200 - F00*F12*W2122 + F02*F12*W2120 + F12*F20*W2102 - F12*F22*W2100 + F00*F11*W2222 - F02*F11*W2220 - F11*F20*W2202 + F11*F22*W2200;
       out_data[5]   =  F21*F21*W1200 - F00*F22*W1121 + F01*F22*W1120 + F20*F22*W1101 - F21*F22*W1100 + F00*F21*W1221 - F01*F21*W1220 - F20*F21*W1201 + F00*F12*W2121 - F01*F12*W2120 - F12*F20*W2101 + F12*F21*W2100 - F00*F11*W2221 + F01*F11*W2220 + F11*F20*W2201 - F11*F21*W2200;
       out_data[6]   =  F01*F22*W1112 - F11*F11*W2202 - F12*F12*W2101 - F02*F22*W1111 - F11*F22*W1102 + F12*F22*W1101 - F01*F21*W1212 + F02*F21*W1211 + F11*F21*W1202 - F12*F21*W1201 - F01*F12*W2112 + F02*F12*W2111 + F11*F12*W2102 + F01*F11*W2212 - F02*F11*W2211 + F11*F12*W2201;
       out_data[7]   =  F12*F12*W2100 - F00*F22*W1112 + F02*F22*W1110 + F10*F22*W1102 - F12*F22*W1100 + F00*F21*W1212 - F02*F21*W1210 - F10*F21*W1202 + F12*F21*W1200 + F00*F12*W2112 - F02*F12*W2110 - F10*F12*W2102 - F00*F11*W2212 + F02*F11*W2210 + F10*F11*W2202 - F11*F12*W2200;
       out_data[8]   =  F11*F11*W2200 + F00*F22*W1111 - F01*F22*W1110 - F10*F22*W1101 + F11*F22*W1100 - F00*F21*W1211 + F01*F21*W1210 + F10*F21*W1201 - F11*F21*W1200 - F00*F12*W2111 + F01*F12*W2110 + F10*F12*W2101 - F11*F12*W2100 + F00*F11*W2211 - F01*F11*W2210 - F10*F11*W2201;
       out_data[9]   =  F12*F22*W1021 - F12*F12*W2021 - F11*F22*W1022 - F22*F22*W1011 + F21*F22*W1012 + F11*F20*W1222 - F12*F20*W1221 - F20*F21*W1212 + F20*F22*W1211 + F11*F12*W2022 - F12*F21*W2012 + F12*F22*W2011 - F10*F11*W2222 + F10*F12*W2221 + F10*F21*W2212 - F10*F22*W2211;
       out_data[10]  =  F22*F22*W1010 + F20*F20*W1212 + F12*F12*W2020 + F10*F10*W2222 + F10*F22*W1022 - F12*F22*W1020 - F20*F22*W1012 - F10*F20*W1222 + F12*F20*W1220 - F20*F22*W1210 - F10*F12*W2022 + F12*F20*W2012 - F12*F22*W2010 - F10*F12*W2220 - F10*F20*W2212 + F10*F22*W2210;
       out_data[11]  =  F11*F22*W1020 - F10*F10*W2221 - F10*F22*W1021 - F20*F20*W1211 + F20*F22*W1011 - F21*F22*W1010 + F10*F20*W1221 - F11*F20*W1220 + F20*F21*W1210 + F10*F12*W2021 - F11*F12*W2020 - F12*F20*W2011 + F12*F21*W2010 + F10*F11*W2220 + F10*F20*W2211 - F10*F21*W2210;
       out_data[12]  =  F22*F22*W1001 + F01*F22*W1022 - F02*F22*W1021 - F21*F22*W1002 - F01*F20*W1222 + F02*F20*W1221 + F20*F21*W1202 - F20*F22*W1201 - F01*F12*W2022 + F02*F12*W2021 + F12*F21*W2002 - F12*F22*W2001 + F01*F10*W2222 - F02*F10*W2221 - F10*F21*W2202 + F10*F22*W2201;
       out_data[13]  =  F02*F22*W1020 - F20*F20*W1202 - F00*F22*W1022 - F22*F22*W1000 + F20*F22*W1002 + F00*F20*W1222 - F02*F20*W1220 + F20*F22*W1200 + F00*F12*W2022 - F02*F12*W2020 - F12*F20*W2002 + F12*F22*W2000 - F00*F10*W2222 + F02*F10*W2220 + F10*F20*W2202 - F10*F22*W2200;
       out_data[14]  =  F20*F20*W1201 + F00*F22*W1021 - F01*F22*W1020 - F20*F22*W1001 + F21*F22*W1000 - F00*F20*W1221 + F01*F20*W1220 - F20*F21*W1200 - F00*F12*W2021 + F01*F12*W2020 + F12*F20*W2001 - F12*F21*W2000 + F00*F10*W2221 - F01*F10*W2220 - F10*F20*W2201 + F10*F21*W2200;
       out_data[15]  =  F12*F12*W2001 - F01*F22*W1012 + F02*F22*W1011 + F11*F22*W1002 - F12*F22*W1001 + F01*F20*W1212 - F02*F20*W1211 - F11*F20*W1202 + F12*F20*W1201 + F01*F12*W2012 - F02*F12*W2011 - F11*F12*W2002 - F01*F10*W2212 + F02*F10*W2211 + F10*F11*W2202 - F10*F12*W2201;
       out_data[16]  =  F00*F22*W1012 - F10*F10*W2202 - F12*F12*W2000 - F02*F22*W1010 - F10*F22*W1002 + F12*F22*W1000 - F00*F20*W1212 + F02*F20*W1210 + F10*F20*W1202 - F12*F20*W1200 - F00*F12*W2012 + F02*F12*W2010 + F10*F12*W2002 + F00*F10*W2212 - F02*F10*W2210 + F10*F12*W2200;
       out_data[17]  =  F10*F10*W2201 - F00*F22*W1011 + F01*F22*W1010 + F10*F22*W1001 - F11*F22*W1000 + F00*F20*W1211 - F01*F20*W1210 - F10*F20*W1201 + F11*F20*W1200 + F00*F12*W2011 - F01*F12*W2010 - F10*F12*W2001 + F11*F12*W2000 - F00*F10*W2211 + F01*F10*W2210 - F10*F11*W2200;
       out_data[18]  =  F11*F21*W1022 - F11*F11*W2022 - F21*F21*W1012 - F12*F21*W1021 + F21*F22*W1011 - F11*F20*W1122 + F12*F20*W1121 + F20*F21*W1112 - F20*F22*W1111 + F11*F12*W2021 + F11*F21*W2012 - F11*F22*W2011 + F10*F11*W2122 - F10*F12*W2121 - F10*F21*W2112 + F10*F22*W2111;
       out_data[19]  =  F12*F21*W1020 - F10*F10*W2122 - F10*F21*W1022 - F20*F20*W1112 + F20*F21*W1012 - F21*F22*W1010 + F10*F20*W1122 - F12*F20*W1120 + F20*F22*W1110 + F10*F11*W2022 - F11*F12*W2020 - F11*F20*W2012 + F11*F22*W2010 + F10*F12*W2120 + F10*F20*W2112 - F10*F22*W2110;
       out_data[20]  =  F21*F21*W1010 + F20*F20*W1111 + F11*F11*W2020 + F10*F10*W2121 + F10*F21*W1021 - F11*F21*W1020 - F20*F21*W1011 - F10*F20*W1121 + F11*F20*W1120 - F20*F21*W1110 - F10*F11*W2021 + F11*F20*W2011 - F11*F21*W2010 - F10*F11*W2120 - F10*F20*W2111 + F10*F21*W2110;
       out_data[21]  =  F21*F21*W1002 - F01*F21*W1022 + F02*F21*W1021 - F21*F22*W1001 + F01*F20*W1122 - F02*F20*W1121 - F20*F21*W1102 + F20*F22*W1101 + F01*F11*W2022 - F02*F11*W2021 - F11*F21*W2002 + F11*F22*W2001 - F01*F10*W2122 + F02*F10*W2121 + F10*F21*W2102 - F10*F22*W2101;
       out_data[22]  =  F20*F20*W1102 + F00*F21*W1022 - F02*F21*W1020 - F20*F21*W1002 + F21*F22*W1000 - F00*F20*W1122 + F02*F20*W1120 - F20*F22*W1100 - F00*F11*W2022 + F02*F11*W2020 + F11*F20*W2002 - F11*F22*W2000 + F00*F10*W2122 - F02*F10*W2120 - F10*F20*W2102 + F10*F22*W2100;
       out_data[23]  =  F01*F21*W1020 - F20*F20*W1101 - F00*F21*W1021 - F21*F21*W1000 + F20*F21*W1001 + F00*F20*W1121 - F01*F20*W1120 + F20*F21*W1100 + F00*F11*W2021 - F01*F11*W2020 - F11*F20*W2001 + F11*F21*W2000 - F00*F10*W2121 + F01*F10*W2120 + F10*F20*W2101 - F10*F21*W2100;
       out_data[24]  =  F11*F11*W2002 + F01*F21*W1012 - F02*F21*W1011 - F11*F21*W1002 + F12*F21*W1001 - F01*F20*W1112 + F02*F20*W1111 + F11*F20*W1102 - F12*F20*W1101 - F01*F11*W2012 + F02*F11*W2011 - F11*F12*W2001 + F01*F10*W2112 - F02*F10*W2111 - F10*F11*W2102 + F10*F12*W2101;
       out_data[25]  =  F10*F10*W2102 - F00*F21*W1012 + F02*F21*W1010 + F10*F21*W1002 - F12*F21*W1000 + F00*F20*W1112 - F02*F20*W1110 - F10*F20*W1102 + F12*F20*W1100 + F00*F11*W2012 - F02*F11*W2010 - F10*F11*W2002 + F11*F12*W2000 - F00*F10*W2112 + F02*F10*W2110 - F10*F12*W2100;
       out_data[26]  =  F00*F21*W1011 - F10*F10*W2101 - F11*F11*W2000 - F01*F21*W1010 - F10*F21*W1001 + F11*F21*W1000 - F00*F20*W1111 + F01*F20*W1110 + F10*F20*W1101 - F11*F20*W1100 - F00*F11*W2011 + F01*F11*W2010 + F10*F11*W2001 + F00*F10*W2111 - F01*F10*W2110 + F10*F11*W2100;
       out_data[27]  =  F12*F22*W0121 - F21*F21*W0212 - F11*F22*W0122 - F22*F22*W0111 + F21*F22*W0112 + F11*F21*W0222 - F12*F21*W0221 + F21*F22*W0211 + F02*F11*W2122 - F02*F12*W2121 - F02*F21*W2112 + F02*F22*W2111 - F01*F11*W2222 + F01*F12*W2221 + F01*F21*W2212 - F01*F22*W2211;
       out_data[28]  =  F22*F22*W0110 + F10*F22*W0122 - F12*F22*W0120 - F20*F22*W0112 - F10*F21*W0222 + F12*F21*W0220 + F20*F21*W0212 - F21*F22*W0210 - F02*F10*W2122 + F02*F12*W2120 + F02*F20*W2112 - F02*F22*W2110 + F01*F10*W2222 - F01*F12*W2220 - F01*F20*W2212 + F01*F22*W2210;
       out_data[29]  =  F21*F21*W0210 - F10*F22*W0121 + F11*F22*W0120 + F20*F22*W0111 - F21*F22*W0110 + F10*F21*W0221 - F11*F21*W0220 - F20*F21*W0211 + F02*F10*W2121 - F02*F11*W2120 - F02*F20*W2111 + F02*F21*W2110 - F01*F10*W2221 + F01*F11*W2220 + F01*F20*W2211 - F01*F21*W2210;
       out_data[30]  =  F22*F22*W0101 + F21*F21*W0202 + F02*F02*W2121 + F01*F01*W2222 + F01*F22*W0122 - F02*F22*W0121 - F21*F22*W0102 - F01*F21*W0222 + F02*F21*W0221 - F21*F22*W0201 - F01*F02*W2122 + F02*F21*W2102 - F02*F22*W2101 - F01*F02*W2221 - F01*F21*W2202 + F01*F22*W2201;
       out_data[31]  =  F02*F22*W0120 - F02*F02*W2120 - F00*F22*W0122 - F22*F22*W0100 + F20*F22*W0102 + F00*F21*W0222 - F02*F21*W0220 - F20*F21*W0202 + F21*F22*W0200 + F00*F02*W2122 - F02*F20*W2102 + F02*F22*W2100 - F00*F01*W2222 + F01*F02*W2220 + F01*F20*W2202 - F01*F22*W2200;
       out_data[32]  =  F00*F22*W0121 - F01*F01*W2220 - F21*F21*W0200 - F01*F22*W0120 - F20*F22*W0101 + F21*F22*W0100 - F00*F21*W0221 + F01*F21*W0220 + F20*F21*W0201 - F00*F02*W2121 + F01*F02*W2120 + F02*F20*W2101 - F02*F21*W2100 + F00*F01*W2221 - F01*F20*W2201 + F01*F21*W2200;
       out_data[33]  =  F02*F22*W0111 - F01*F01*W2212 - F01*F22*W0112 - F02*F02*W2111 + F11*F22*W0102 - F12*F22*W0101 + F01*F21*W0212 - F02*F21*W0211 - F11*F21*W0202 + F12*F21*W0201 + F01*F02*W2112 - F02*F11*W2102 + F02*F12*W2101 + F01*F02*W2211 + F01*F11*W2202 - F01*F12*W2201;
       out_data[34]  =  F02*F02*W2110 + F00*F22*W0112 - F02*F22*W0110 - F10*F22*W0102 + F12*F22*W0100 - F00*F21*W0212 + F02*F21*W0210 + F10*F21*W0202 - F12*F21*W0200 - F00*F02*W2112 + F02*F10*W2102 - F02*F12*W2100 + F00*F01*W2212 - F01*F02*W2210 - F01*F10*W2202 + F01*F12*W2200;
       out_data[35]  =  F01*F01*W2210 - F00*F22*W0111 + F01*F22*W0110 + F10*F22*W0101 - F11*F22*W0100 + F00*F21*W0211 - F01*F21*W0210 - F10*F21*W0201 + F11*F21*W0200 + F00*F02*W2111 - F01*F02*W2110 - F02*F10*W2101 + F02*F11*W2100 - F00*F01*W2211 + F01*F10*W2201 - F01*F11*W2200;
       out_data[36]  =  F22*F22*W0011 + F11*F22*W0022 - F12*F22*W0021 - F21*F22*W0012 - F11*F20*W0222 + F12*F20*W0221 + F20*F21*W0212 - F20*F22*W0211 - F02*F11*W2022 + F02*F12*W2021 + F02*F21*W2012 - F02*F22*W2011 + F00*F11*W2222 - F00*F12*W2221 - F00*F21*W2212 + F00*F22*W2211;
       out_data[37]  =  F12*F22*W0020 - F20*F20*W0212 - F10*F22*W0022 - F22*F22*W0010 + F20*F22*W0012 + F10*F20*W0222 - F12*F20*W0220 + F20*F22*W0210 + F02*F10*W2022 - F02*F12*W2020 - F02*F20*W2012 + F02*F22*W2010 - F00*F10*W2222 + F00*F12*W2220 + F00*F20*W2212 - F00*F22*W2210;
       out_data[38]  =  F20*F20*W0211 + F10*F22*W0021 - F11*F22*W0020 - F20*F22*W0011 + F21*F22*W0010 - F10*F20*W0221 + F11*F20*W0220 - F20*F21*W0210 - F02*F10*W2021 + F02*F11*W2020 + F02*F20*W2011 - F02*F21*W2010 + F00*F10*W2221 - F00*F11*W2220 - F00*F20*W2211 + F00*F21*W2210;
       out_data[39]  =  F02*F22*W0021 - F02*F02*W2021 - F01*F22*W0022 - F22*F22*W0001 + F21*F22*W0002 + F01*F20*W0222 - F02*F20*W0221 - F20*F21*W0202 + F20*F22*W0201 + F01*F02*W2022 - F02*F21*W2002 + F02*F22*W2001 - F00*F01*W2222 + F00*F02*W2221 + F00*F21*W2202 - F00*F22*W2201;
       out_data[40]  =  F22*F22*W0000 + F20*F20*W0202 + F02*F02*W2020 + F00*F00*W2222 + F00*F22*W0022 - F02*F22*W0020 - F20*F22*W0002 - F00*F20*W0222 + F02*F20*W0220 - F20*F22*W0200 - F00*F02*W2022 + F02*F20*W2002 - F02*F22*W2000 - F00*F02*W2220 - F00*F20*W2202 + F00*F22*W2200;
       out_data[41]  =  F01*F22*W0020 - F00*F00*W2221 - F00*F22*W0021 - F20*F20*W0201 + F20*F22*W0001 - F21*F22*W0000 + F00*F20*W0221 - F01*F20*W0220 + F20*F21*W0200 + F00*F02*W2021 - F01*F02*W2020 - F02*F20*W2001 + F02*F21*W2000 + F00*F01*W2220 + F00*F20*W2201 - F00*F21*W2200;
       out_data[42]  =  F02*F02*W2011 + F01*F22*W0012 - F02*F22*W0011 - F11*F22*W0002 + F12*F22*W0001 - F01*F20*W0212 + F02*F20*W0211 + F11*F20*W0202 - F12*F20*W0201 - F01*F02*W2012 + F02*F11*W2002 - F02*F12*W2001 + F00*F01*W2212 - F00*F02*W2211 - F00*F11*W2202 + F00*F12*W2201;
       out_data[43]  =  F02*F22*W0010 - F00*F00*W2212 - F00*F22*W0012 - F02*F02*W2010 + F10*F22*W0002 - F12*F22*W0000 + F00*F20*W0212 - F02*F20*W0210 - F10*F20*W0202 + F12*F20*W0200 + F00*F02*W2012 - F02*F10*W2002 + F02*F12*W2000 + F00*F02*W2210 + F00*F10*W2202 - F00*F12*W2200;
       out_data[44]  =  F00*F00*W2211 + F00*F22*W0011 - F01*F22*W0010 - F10*F22*W0001 + F11*F22*W0000 - F00*F20*W0211 + F01*F20*W0210 + F10*F20*W0201 - F11*F20*W0200 - F00*F02*W2011 + F01*F02*W2010 + F02*F10*W2001 - F02*F11*W2000 - F00*F01*W2210 - F00*F10*W2201 + F00*F11*W2200;
       out_data[45]  =  F21*F21*W0012 - F11*F21*W0022 + F12*F21*W0021 - F21*F22*W0011 + F11*F20*W0122 - F12*F20*W0121 - F20*F21*W0112 + F20*F22*W0111 + F01*F11*W2022 - F01*F12*W2021 - F01*F21*W2012 + F01*F22*W2011 - F00*F11*W2122 + F00*F12*W2121 + F00*F21*W2112 - F00*F22*W2111;
       out_data[46]  =  F20*F20*W0112 + F10*F21*W0022 - F12*F21*W0020 - F20*F21*W0012 + F21*F22*W0010 - F10*F20*W0122 + F12*F20*W0120 - F20*F22*W0110 - F01*F10*W2022 + F01*F12*W2020 + F01*F20*W2012 - F01*F22*W2010 + F00*F10*W2122 - F00*F12*W2120 - F00*F20*W2112 + F00*F22*W2110;
       out_data[47]  =  F11*F21*W0020 - F20*F20*W0111 - F10*F21*W0021 - F21*F21*W0010 + F20*F21*W0011 + F10*F20*W0121 - F11*F20*W0120 + F20*F21*W0110 + F01*F10*W2021 - F01*F11*W2020 - F01*F20*W2011 + F01*F21*W2010 - F00*F10*W2121 + F00*F11*W2120 + F00*F20*W2111 - F00*F21*W2110;
       out_data[48]  =  F01*F21*W0022 - F01*F01*W2022 - F21*F21*W0002 - F02*F21*W0021 + F21*F22*W0001 - F01*F20*W0122 + F02*F20*W0121 + F20*F21*W0102 - F20*F22*W0101 + F01*F02*W2021 + F01*F21*W2002 - F01*F22*W2001 + F00*F01*W2122 - F00*F02*W2121 - F00*F21*W2102 + F00*F22*W2101;
       out_data[49]  =  F02*F21*W0020 - F00*F00*W2122 - F00*F21*W0022 - F20*F20*W0102 + F20*F21*W0002 - F21*F22*W0000 + F00*F20*W0122 - F02*F20*W0120 + F20*F22*W0100 + F00*F01*W2022 - F01*F02*W2020 - F01*F20*W2002 + F01*F22*W2000 + F00*F02*W2120 + F00*F20*W2102 - F00*F22*W2100;
       out_data[50]  =  F21*F21*W0000 + F20*F20*W0101 + F01*F01*W2020 + F00*F00*W2121 + F00*F21*W0021 - F01*F21*W0020 - F20*F21*W0001 - F00*F20*W0121 + F01*F20*W0120 - F20*F21*W0100 - F00*F01*W2021 + F01*F20*W2001 - F01*F21*W2000 - F00*F01*W2120 - F00*F20*W2101 + F00*F21*W2100;
       out_data[51]  =  F01*F01*W2012 - F01*F21*W0012 + F02*F21*W0011 + F11*F21*W0002 - F12*F21*W0001 + F01*F20*W0112 - F02*F20*W0111 - F11*F20*W0102 + F12*F20*W0101 - F01*F02*W2011 - F01*F11*W2002 + F01*F12*W2001 - F00*F01*W2112 + F00*F02*W2111 + F00*F11*W2102 - F00*F12*W2101;
       out_data[52]  =  F00*F00*W2112 + F00*F21*W0012 - F02*F21*W0010 - F10*F21*W0002 + F12*F21*W0000 - F00*F20*W0112 + F02*F20*W0110 + F10*F20*W0102 - F12*F20*W0100 - F00*F01*W2012 + F01*F02*W2010 + F01*F10*W2002 - F01*F12*W2000 - F00*F02*W2110 - F00*F10*W2102 + F00*F12*W2100;
       out_data[53]  =  F01*F21*W0010 - F00*F00*W2111 - F00*F21*W0011 - F01*F01*W2010 + F10*F21*W0001 - F11*F21*W0000 + F00*F20*W0111 - F01*F20*W0110 - F10*F20*W0101 + F11*F20*W0100 + F00*F01*W2011 - F01*F10*W2001 + F01*F11*W2000 + F00*F01*W2110 + F00*F10*W2101 - F00*F11*W2100;
       out_data[54]  =  F11*F12*W0122 - F11*F11*W0222 - F12*F12*W0121 - F12*F21*W0112 + F12*F22*W0111 + F11*F12*W0221 + F11*F21*W0212 - F11*F22*W0211 - F02*F11*W1122 + F02*F12*W1121 + F02*F21*W1112 - F02*F22*W1111 + F01*F11*W1222 - F01*F12*W1221 - F01*F21*W1212 + F01*F22*W1211;
       out_data[55]  =  F12*F12*W0120 - F10*F12*W0122 + F12*F20*W0112 - F12*F22*W0110 + F10*F11*W0222 - F11*F12*W0220 - F11*F20*W0212 + F11*F22*W0210 + F02*F10*W1122 - F02*F12*W1120 - F02*F20*W1112 + F02*F22*W1110 - F01*F10*W1222 + F01*F12*W1220 + F01*F20*W1212 - F01*F22*W1210;
       out_data[56]  =  F11*F11*W0220 + F10*F12*W0121 - F11*F12*W0120 - F12*F20*W0111 + F12*F21*W0110 - F10*F11*W0221 + F11*F20*W0211 - F11*F21*W0210 - F02*F10*W1121 + F02*F11*W1120 + F02*F20*W1111 - F02*F21*W1110 + F01*F10*W1221 - F01*F11*W1220 - F01*F20*W1211 + F01*F21*W1210;
       out_data[57]  =  F02*F12*W0121 - F01*F01*W1222 - F01*F12*W0122 - F02*F02*W1121 + F12*F21*W0102 - F12*F22*W0101 + F01*F11*W0222 - F02*F11*W0221 - F11*F21*W0202 + F11*F22*W0201 + F01*F02*W1122 - F02*F21*W1102 + F02*F22*W1101 + F01*F02*W1221 + F01*F21*W1202 - F01*F22*W1201;
       out_data[58]  =  F02*F02*W1120 + F00*F12*W0122 - F02*F12*W0120 - F12*F20*W0102 + F12*F22*W0100 - F00*F11*W0222 + F02*F11*W0220 + F11*F20*W0202 - F11*F22*W0200 - F00*F02*W1122 + F02*F20*W1102 - F02*F22*W1100 + F00*F01*W1222 - F01*F02*W1220 - F01*F20*W1202 + F01*F22*W1200;
       out_data[59]  =  F01*F01*W1220 - F00*F12*W0121 + F01*F12*W0120 + F12*F20*W0101 - F12*F21*W0100 + F00*F11*W0221 - F01*F11*W0220 - F11*F20*W0201 + F11*F21*W0200 + F00*F02*W1121 - F01*F02*W1120 - F02*F20*W1101 + F02*F21*W1100 - F00*F01*W1221 + F01*F20*W1201 - F01*F21*W1200;
       out_data[60]  =  F12*F12*W0101 + F11*F11*W0202 + F02*F02*W1111 + F01*F01*W1212 + F01*F12*W0112 - F02*F12*W0111 - F11*F12*W0102 - F01*F11*W0212 + F02*F11*W0211 - F11*F12*W0201 - F01*F02*W1112 + F02*F11*W1102 - F02*F12*W1101 - F01*F02*W1211 - F01*F11*W1202 + F01*F12*W1201;
       out_data[61]  =  F02*F12*W0110 - F02*F02*W1110 - F00*F12*W0112 - F12*F12*W0100 + F10*F12*W0102 + F00*F11*W0212 - F02*F11*W0210 - F10*F11*W0202 + F11*F12*W0200 + F00*F02*W1112 - F02*F10*W1102 + F02*F12*W1100 - F00*F01*W1212 + F01*F02*W1210 + F01*F10*W1202 - F01*F12*W1200;
       out_data[62]  =  F00*F12*W0111 - F01*F01*W1210 - F11*F11*W0200 - F01*F12*W0110 - F10*F12*W0101 + F11*F12*W0100 - F00*F11*W0211 + F01*F11*W0210 + F10*F11*W0201 - F00*F02*W1111 + F01*F02*W1110 + F02*F10*W1101 - F02*F11*W1100 + F00*F01*W1211 - F01*F10*W1201 + F01*F11*W1200;
       out_data[63]  =  F12*F12*W0021 - F11*F12*W0022 + F12*F21*W0012 - F12*F22*W0011 + F10*F11*W0222 - F10*F12*W0221 - F10*F21*W0212 + F10*F22*W0211 + F02*F11*W1022 - F02*F12*W1021 - F02*F21*W1012 + F02*F22*W1011 - F00*F11*W1222 + F00*F12*W1221 + F00*F21*W1212 - F00*F22*W1211;
       out_data[64]  =  F10*F12*W0022 - F10*F10*W0222 - F12*F12*W0020 - F12*F20*W0012 + F12*F22*W0010 + F10*F12*W0220 + F10*F20*W0212 - F10*F22*W0210 - F02*F10*W1022 + F02*F12*W1020 + F02*F20*W1012 - F02*F22*W1010 + F00*F10*W1222 - F00*F12*W1220 - F00*F20*W1212 + F00*F22*W1210;
       out_data[65]  =  F10*F10*W0221 - F10*F12*W0021 + F11*F12*W0020 + F12*F20*W0011 - F12*F21*W0010 - F10*F11*W0220 - F10*F20*W0211 + F10*F21*W0210 + F02*F10*W1021 - F02*F11*W1020 - F02*F20*W1011 + F02*F21*W1010 - F00*F10*W1221 + F00*F11*W1220 + F00*F20*W1211 - F00*F21*W1210;
       out_data[66]  =  F02*F02*W1021 + F01*F12*W0022 - F02*F12*W0021 - F12*F21*W0002 + F12*F22*W0001 - F01*F10*W0222 + F02*F10*W0221 + F10*F21*W0202 - F10*F22*W0201 - F01*F02*W1022 + F02*F21*W1002 - F02*F22*W1001 + F00*F01*W1222 - F00*F02*W1221 - F00*F21*W1202 + F00*F22*W1201;
       out_data[67]  =  F02*F12*W0020 - F00*F00*W1222 - F00*F12*W0022 - F02*F02*W1020 + F12*F20*W0002 - F12*F22*W0000 + F00*F10*W0222 - F02*F10*W0220 - F10*F20*W0202 + F10*F22*W0200 + F00*F02*W1022 - F02*F20*W1002 + F02*F22*W1000 + F00*F02*W1220 + F00*F20*W1202 - F00*F22*W1200;
       out_data[68]  =  F00*F00*W1221 + F00*F12*W0021 - F01*F12*W0020 - F12*F20*W0001 + F12*F21*W0000 - F00*F10*W0221 + F01*F10*W0220 + F10*F20*W0201 - F10*F21*W0200 - F00*F02*W1021 + F01*F02*W1020 + F02*F20*W1001 - F02*F21*W1000 - F00*F01*W1220 - F00*F20*W1201 + F00*F21*W1200;
       out_data[69]  =  F02*F12*W0011 - F02*F02*W1011 - F01*F12*W0012 - F12*F12*W0001 + F11*F12*W0002 + F01*F10*W0212 - F02*F10*W0211 - F10*F11*W0202 + F10*F12*W0201 + F01*F02*W1012 - F02*F11*W1002 + F02*F12*W1001 - F00*F01*W1212 + F00*F02*W1211 + F00*F11*W1202 - F00*F12*W1201;
       out_data[70]  =  F12*F12*W0000 + F10*F10*W0202 + F02*F02*W1010 + F00*F00*W1212 + F00*F12*W0012 - F02*F12*W0010 - F10*F12*W0002 - F00*F10*W0212 + F02*F10*W0210 - F10*F12*W0200 - F00*F02*W1012 + F02*F10*W1002 - F02*F12*W1000 - F00*F02*W1210 - F00*F10*W1202 + F00*F12*W1200;
       out_data[71]  =  F01*F12*W0010 - F00*F00*W1211 - F00*F12*W0011 - F10*F10*W0201 + F10*F12*W0001 - F11*F12*W0000 + F00*F10*W0211 - F01*F10*W0210 + F10*F11*W0200 + F00*F02*W1011 - F01*F02*W1010 - F02*F10*W1001 + F02*F11*W1000 + F00*F01*W1210 + F00*F10*W1201 - F00*F11*W1200;
       out_data[72]  =  F11*F11*W0022 - F11*F12*W0021 - F11*F21*W0012 + F11*F22*W0011 - F10*F11*W0122 + F10*F12*W0121 + F10*F21*W0112 - F10*F22*W0111 - F01*F11*W1022 + F01*F12*W1021 + F01*F21*W1012 - F01*F22*W1011 + F00*F11*W1122 - F00*F12*W1121 - F00*F21*W1112 + F00*F22*W1111;
       out_data[73]  =  F10*F10*W0122 - F10*F11*W0022 + F11*F12*W0020 + F11*F20*W0012 - F11*F22*W0010 - F10*F12*W0120 - F10*F20*W0112 + F10*F22*W0110 + F01*F10*W1022 - F01*F12*W1020 - F01*F20*W1012 + F01*F22*W1010 - F00*F10*W1122 + F00*F12*W1120 + F00*F20*W1112 - F00*F22*W1110;
       out_data[74]  =  F10*F11*W0021 - F10*F10*W0121 - F11*F11*W0020 - F11*F20*W0011 + F11*F21*W0010 + F10*F11*W0120 + F10*F20*W0111 - F10*F21*W0110 - F01*F10*W1021 + F01*F11*W1020 + F01*F20*W1011 - F01*F21*W1010 + F00*F10*W1121 - F00*F11*W1120 - F00*F20*W1111 + F00*F21*W1110;
       out_data[75]  =  F01*F01*W1022 - F01*F11*W0022 + F02*F11*W0021 + F11*F21*W0002 - F11*F22*W0001 + F01*F10*W0122 - F02*F10*W0121 - F10*F21*W0102 + F10*F22*W0101 - F01*F02*W1021 - F01*F21*W1002 + F01*F22*W1001 - F00*F01*W1122 + F00*F02*W1121 + F00*F21*W1102 - F00*F22*W1101;
       out_data[76]  =  F00*F00*W1122 + F00*F11*W0022 - F02*F11*W0020 - F11*F20*W0002 + F11*F22*W0000 - F00*F10*W0122 + F02*F10*W0120 + F10*F20*W0102 - F10*F22*W0100 - F00*F01*W1022 + F01*F02*W1020 + F01*F20*W1002 - F01*F22*W1000 - F00*F02*W1120 - F00*F20*W1102 + F00*F22*W1100;
       out_data[77]  =  F01*F11*W0020 - F00*F00*W1121 - F00*F11*W0021 - F01*F01*W1020 + F11*F20*W0001 - F11*F21*W0000 + F00*F10*W0121 - F01*F10*W0120 - F10*F20*W0101 + F10*F21*W0100 + F00*F01*W1021 - F01*F20*W1001 + F01*F21*W1000 + F00*F01*W1120 + F00*F20*W1101 - F00*F21*W1100;
       out_data[78]  =  F01*F11*W0012 - F01*F01*W1012 - F11*F11*W0002 - F02*F11*W0011 + F11*F12*W0001 - F01*F10*W0112 + F02*F10*W0111 + F10*F11*W0102 - F10*F12*W0101 + F01*F02*W1011 + F01*F11*W1002 - F01*F12*W1001 + F00*F01*W1112 - F00*F02*W1111 - F00*F11*W1102 + F00*F12*W1101;
       out_data[79]  =  F02*F11*W0010 - F00*F00*W1112 - F00*F11*W0012 - F10*F10*W0102 + F10*F11*W0002 - F11*F12*W0000 + F00*F10*W0112 - F02*F10*W0110 + F10*F12*W0100 + F00*F01*W1012 - F01*F02*W1010 - F01*F10*W1002 + F01*F12*W1000 + F00*F02*W1110 + F00*F10*W1102 - F00*F12*W1100;
       out_data[80]  =  F11*F11*W0000 + F10*F10*W0101 + F01*F01*W1010 + F00*F00*W1111 + F00*F11*W0011 - F01*F11*W0010 - F10*F11*W0001 - F00*F10*W0111 + F01*F10*W0110 - F10*F11*W0100 - F00*F01*W1011 + F01*F10*W1001 - F01*F11*W1000 - F00*F01*W1110 - F00*F10*W1101 + F00*F11*W1100;
    }
    else{
       T W0000  =  H_data[0];   T W0001  =  H_data[1];   T W0010  =  H_data[2];   T W0011  =  H_data[3];   T W0100  =  H_data[4];
       T W0101  =  H_data[5];   T W0110  =  H_data[6];   T W0111  =  H_data[7];   T W1000  =  H_data[8];   T W1001  =  H_data[9];
       T W1010  =  H_data[10];  T W1011  =  H_data[11];  T W1100  =  H_data[12];  T W1101  =  H_data[13];  T W1110  =  H_data[14];
       T W1111  =  H_data[15];


       out_data[0]   =   W1111;   out_data[1]   =  -W1110;    out_data[2]   =  -W1101;   out_data[3]   =   W1100;
       out_data[4]   =  -W1011;   out_data[5]   =   W1010;    out_data[6]   =   W1001;   out_data[7]   =  -W1000;
       out_data[8]   =  -W0111;   out_data[9]   =   W0110;    out_data[10]  =   W0101;   out_data[11]  =  -W0100;
       out_data[12]  =   W0011;   out_data[13]  =  -W0010;    out_data[14]  =  -W0001;   out_data[15]  =   W0000;
    }
}

tensor<real,5> H_HH(const tensor<real,3> &F, const tensor<real,5> &WHH,  const int &ndim) {
    //--------------------------------------------------------------
    // Obtain FH contribution in the Hessian for every Gauss point
    //--------------------------------------------------------------
    auto ngauss        =  F.dimension(0);
    tensor<real,5> out(ngauss,ndim,ndim,ndim,ndim);
    for (auto i=0; i<ngauss; ++i) {
        tensor<real,2> F_2d   =  F.chip(i,0);
        tensor<real,4> WHH_4d =  WHH.chip(i,0);
        __Hessian_HH_contribution__(F_2d.data(),WHH_4d.data(), out.data()+ndim*ndim*ndim*ndim*i,ndim);
    }
    return out;
}


//template<typename T>
//SC_INLINE void __Hessian_HJ_contribution__(const T *__restrict__ F_data, const T *__restrict__ H_data, const T *__restrict__ WHJ_data, T *__restrict__ out_data, const int &ndim) {
//    //--------------------------------------------------------------
//    // Contribution HF in the Hessian operator
//    //--------------------------------------------------------------
//    if (ndim==3){
//       T F00         =  F_data[0];
//       T F11         =  F_data[4];
//       T F22         =  F_data[8];
//       T F01         =  F_data[1];
//       T F02         =  F_data[2];
//       T F12         =  F_data[5];
//       T F10         =  F_data[3];
//       T F20         =  F_data[6];
//       T F21         =  F_data[7];

//       T H00         =  H_data[0];
//       T H11         =  H_data[4];
//       T H22         =  H_data[8];
//       T H01         =  H_data[1];
//       T H02         =  H_data[2];
//       T H12         =  H_data[5];
//       T H10         =  H_data[3];
//       T H20         =  H_data[6];
//       T H21         =  H_data[7];

//       T WHJ00       =  WHJ_data[0];
//       T WHJ11       =  WHJ_data[4];
//       T WHJ22       =  WHJ_data[8];
//       T WHJ01       =  WHJ_data[1];
//       T WHJ02       =  WHJ_data[2];
//       T WHJ12       =  WHJ_data[5];
//       T WHJ10       =  WHJ_data[3];
//       T WHJ20       =  WHJ_data[6];
//       T WHJ21       =  WHJ_data[7];

//       out_data[0]   =  F11*H00*WHJ22 - F12*H00*WHJ21 - F21*H00*WHJ12 + F22*H00*WHJ11;
//       out_data[1]   =  F11*H01*WHJ22 - F12*H01*WHJ21 - F21*H01*WHJ12 + F22*H01*WHJ11;
//       out_data[2]   =  F11*H02*WHJ22 - F12*H02*WHJ21 - F21*H02*WHJ12 + F22*H02*WHJ11;
//       out_data[3]   =  F11*H10*WHJ22 - F12*H10*WHJ21 - F21*H10*WHJ12 + F22*H10*WHJ11;
//       out_data[4]   =  F11*H11*WHJ22 - F12*H11*WHJ21 - F21*H11*WHJ12 + F22*H11*WHJ11;
//       out_data[5]   =  F11*H12*WHJ22 - F12*H12*WHJ21 - F21*H12*WHJ12 + F22*H12*WHJ11;
//       out_data[6]   =  F11*H20*WHJ22 - F12*H20*WHJ21 - F21*H20*WHJ12 + F22*H20*WHJ11;
//       out_data[7]   =  F11*H21*WHJ22 - F12*H21*WHJ21 - F21*H21*WHJ12 + F22*H21*WHJ11;
//       out_data[8]   =  F11*H22*WHJ22 - F12*H22*WHJ21 - F21*H22*WHJ12 + F22*H22*WHJ11;
//       out_data[9]   =  F12*H00*WHJ20 - F10*H00*WHJ22 + F20*H00*WHJ12 - F22*H00*WHJ10;
//       out_data[10]  =  F12*H01*WHJ20 - F10*H01*WHJ22 + F20*H01*WHJ12 - F22*H01*WHJ10;
//       out_data[11]  =  F12*H02*WHJ20 - F10*H02*WHJ22 + F20*H02*WHJ12 - F22*H02*WHJ10;
//       out_data[12]  =  F12*H10*WHJ20 - F10*H10*WHJ22 + F20*H10*WHJ12 - F22*H10*WHJ10;
//       out_data[13]  =  F12*H11*WHJ20 - F10*H11*WHJ22 + F20*H11*WHJ12 - F22*H11*WHJ10;
//       out_data[14]  =  F12*H12*WHJ20 - F10*H12*WHJ22 + F20*H12*WHJ12 - F22*H12*WHJ10;
//       out_data[15]  =  F12*H20*WHJ20 - F10*H20*WHJ22 + F20*H20*WHJ12 - F22*H20*WHJ10;
//       out_data[16]  =  F12*H21*WHJ20 - F10*H21*WHJ22 + F20*H21*WHJ12 - F22*H21*WHJ10;
//       out_data[17]  =  F12*H22*WHJ20 - F10*H22*WHJ22 + F20*H22*WHJ12 - F22*H22*WHJ10;
//       out_data[18]  =  F10*H00*WHJ21 - F11*H00*WHJ20 - F20*H00*WHJ11 + F21*H00*WHJ10;
//       out_data[19]  =  F10*H01*WHJ21 - F11*H01*WHJ20 - F20*H01*WHJ11 + F21*H01*WHJ10;
//       out_data[20]  =  F10*H02*WHJ21 - F11*H02*WHJ20 - F20*H02*WHJ11 + F21*H02*WHJ10;
//       out_data[21]  =  F10*H10*WHJ21 - F11*H10*WHJ20 - F20*H10*WHJ11 + F21*H10*WHJ10;
//       out_data[22]  =  F10*H11*WHJ21 - F11*H11*WHJ20 - F20*H11*WHJ11 + F21*H11*WHJ10;
//       out_data[23]  =  F10*H12*WHJ21 - F11*H12*WHJ20 - F20*H12*WHJ11 + F21*H12*WHJ10;
//       out_data[24]  =  F10*H20*WHJ21 - F11*H20*WHJ20 - F20*H20*WHJ11 + F21*H20*WHJ10;
//       out_data[25]  =  F10*H21*WHJ21 - F11*H21*WHJ20 - F20*H21*WHJ11 + F21*H21*WHJ10;
//       out_data[26]  =  F10*H22*WHJ21 - F11*H22*WHJ20 - F20*H22*WHJ11 + F21*H22*WHJ10;
//       out_data[27]  =  F02*H00*WHJ21 - F01*H00*WHJ22 + F21*H00*WHJ02 - F22*H00*WHJ01;
//       out_data[28]  =  F02*H01*WHJ21 - F01*H01*WHJ22 + F21*H01*WHJ02 - F22*H01*WHJ01;
//       out_data[29]  =  F02*H02*WHJ21 - F01*H02*WHJ22 + F21*H02*WHJ02 - F22*H02*WHJ01;
//       out_data[30]  =  F02*H10*WHJ21 - F01*H10*WHJ22 + F21*H10*WHJ02 - F22*H10*WHJ01;
//       out_data[31]  =  F02*H11*WHJ21 - F01*H11*WHJ22 + F21*H11*WHJ02 - F22*H11*WHJ01;
//       out_data[32]  =  F02*H12*WHJ21 - F01*H12*WHJ22 + F21*H12*WHJ02 - F22*H12*WHJ01;
//       out_data[33]  =  F02*H20*WHJ21 - F01*H20*WHJ22 + F21*H20*WHJ02 - F22*H20*WHJ01;
//       out_data[34]  =  F02*H21*WHJ21 - F01*H21*WHJ22 + F21*H21*WHJ02 - F22*H21*WHJ01;
//       out_data[35]  =  F02*H22*WHJ21 - F01*H22*WHJ22 + F21*H22*WHJ02 - F22*H22*WHJ01;
//       out_data[36]  =  F00*H00*WHJ22 - F02*H00*WHJ20 - F20*H00*WHJ02 + F22*H00*WHJ00;
//       out_data[37]  =  F00*H01*WHJ22 - F02*H01*WHJ20 - F20*H01*WHJ02 + F22*H01*WHJ00;
//       out_data[38]  =  F00*H02*WHJ22 - F02*H02*WHJ20 - F20*H02*WHJ02 + F22*H02*WHJ00;
//       out_data[39]  =  F00*H10*WHJ22 - F02*H10*WHJ20 - F20*H10*WHJ02 + F22*H10*WHJ00;
//       out_data[40]  =  F00*H11*WHJ22 - F02*H11*WHJ20 - F20*H11*WHJ02 + F22*H11*WHJ00;
//       out_data[41]  =  F00*H12*WHJ22 - F02*H12*WHJ20 - F20*H12*WHJ02 + F22*H12*WHJ00;
//       out_data[42]  =  F00*H20*WHJ22 - F02*H20*WHJ20 - F20*H20*WHJ02 + F22*H20*WHJ00;
//       out_data[43]  =  F00*H21*WHJ22 - F02*H21*WHJ20 - F20*H21*WHJ02 + F22*H21*WHJ00;
//       out_data[44]  =  F00*H22*WHJ22 - F02*H22*WHJ20 - F20*H22*WHJ02 + F22*H22*WHJ00;
//       out_data[45]  =  F01*H00*WHJ20 - F00*H00*WHJ21 + F20*H00*WHJ01 - F21*H00*WHJ00;
//       out_data[46]  =  F01*H01*WHJ20 - F00*H01*WHJ21 + F20*H01*WHJ01 - F21*H01*WHJ00;
//       out_data[47]  =  F01*H02*WHJ20 - F00*H02*WHJ21 + F20*H02*WHJ01 - F21*H02*WHJ00;
//       out_data[48]  =  F01*H10*WHJ20 - F00*H10*WHJ21 + F20*H10*WHJ01 - F21*H10*WHJ00;
//       out_data[49]  =  F01*H11*WHJ20 - F00*H11*WHJ21 + F20*H11*WHJ01 - F21*H11*WHJ00;
//       out_data[50]  =  F01*H12*WHJ20 - F00*H12*WHJ21 + F20*H12*WHJ01 - F21*H12*WHJ00;
//       out_data[51]  =  F01*H20*WHJ20 - F00*H20*WHJ21 + F20*H20*WHJ01 - F21*H20*WHJ00;
//       out_data[52]  =  F01*H21*WHJ20 - F00*H21*WHJ21 + F20*H21*WHJ01 - F21*H21*WHJ00;
//       out_data[53]  =  F01*H22*WHJ20 - F00*H22*WHJ21 + F20*H22*WHJ01 - F21*H22*WHJ00;
//       out_data[54]  =  F01*H00*WHJ12 - F02*H00*WHJ11 - F11*H00*WHJ02 + F12*H00*WHJ01;
//       out_data[55]  =  F01*H01*WHJ12 - F02*H01*WHJ11 - F11*H01*WHJ02 + F12*H01*WHJ01;
//       out_data[56]  =  F01*H02*WHJ12 - F02*H02*WHJ11 - F11*H02*WHJ02 + F12*H02*WHJ01;
//       out_data[57]  =  F01*H10*WHJ12 - F02*H10*WHJ11 - F11*H10*WHJ02 + F12*H10*WHJ01;
//       out_data[58]  =  F01*H11*WHJ12 - F02*H11*WHJ11 - F11*H11*WHJ02 + F12*H11*WHJ01;
//       out_data[59]  =  F01*H12*WHJ12 - F02*H12*WHJ11 - F11*H12*WHJ02 + F12*H12*WHJ01;
//       out_data[60]  =  F01*H20*WHJ12 - F02*H20*WHJ11 - F11*H20*WHJ02 + F12*H20*WHJ01;
//       out_data[61]  =  F01*H21*WHJ12 - F02*H21*WHJ11 - F11*H21*WHJ02 + F12*H21*WHJ01;
//       out_data[62]  =  F01*H22*WHJ12 - F02*H22*WHJ11 - F11*H22*WHJ02 + F12*H22*WHJ01;
//       out_data[63]  =  F02*H00*WHJ10 - F00*H00*WHJ12 + F10*H00*WHJ02 - F12*H00*WHJ00;
//       out_data[64]  =  F02*H01*WHJ10 - F00*H01*WHJ12 + F10*H01*WHJ02 - F12*H01*WHJ00;
//       out_data[65]  =  F02*H02*WHJ10 - F00*H02*WHJ12 + F10*H02*WHJ02 - F12*H02*WHJ00;
//       out_data[66]  =  F02*H10*WHJ10 - F00*H10*WHJ12 + F10*H10*WHJ02 - F12*H10*WHJ00;
//       out_data[67]  =  F02*H11*WHJ10 - F00*H11*WHJ12 + F10*H11*WHJ02 - F12*H11*WHJ00;
//       out_data[68]  =  F02*H12*WHJ10 - F00*H12*WHJ12 + F10*H12*WHJ02 - F12*H12*WHJ00;
//       out_data[69]  =  F02*H20*WHJ10 - F00*H20*WHJ12 + F10*H20*WHJ02 - F12*H20*WHJ00;
//       out_data[70]  =  F02*H21*WHJ10 - F00*H21*WHJ12 + F10*H21*WHJ02 - F12*H21*WHJ00;
//       out_data[71]  =  F02*H22*WHJ10 - F00*H22*WHJ12 + F10*H22*WHJ02 - F12*H22*WHJ00;
//       out_data[72]  =  F00*H00*WHJ11 - F01*H00*WHJ10 - F10*H00*WHJ01 + F11*H00*WHJ00;
//       out_data[73]  =  F00*H01*WHJ11 - F01*H01*WHJ10 - F10*H01*WHJ01 + F11*H01*WHJ00;
//       out_data[74]  =  F00*H02*WHJ11 - F01*H02*WHJ10 - F10*H02*WHJ01 + F11*H02*WHJ00;
//       out_data[75]  =  F00*H10*WHJ11 - F01*H10*WHJ10 - F10*H10*WHJ01 + F11*H10*WHJ00;
//       out_data[76]  =  F00*H11*WHJ11 - F01*H11*WHJ10 - F10*H11*WHJ01 + F11*H11*WHJ00;
//       out_data[77]  =  F00*H12*WHJ11 - F01*H12*WHJ10 - F10*H12*WHJ01 + F11*H12*WHJ00;
//       out_data[78]  =  F00*H20*WHJ11 - F01*H20*WHJ10 - F10*H20*WHJ01 + F11*H20*WHJ00;
//       out_data[79]  =  F00*H21*WHJ11 - F01*H21*WHJ10 - F10*H21*WHJ01 + F11*H21*WHJ00;
//       out_data[80]  =  F00*H22*WHJ11 - F01*H22*WHJ10 - F10*H22*WHJ01 + F11*H22*WHJ00;
//    }
//    else{
//       T H00         =  H_data[0];
//       T H11         =  H_data[3];
//       T H01         =  H_data[1];
//       T H10         =  H_data[2];

//       T WHJ00       =  WHJ_data[0];
//       T WHJ11       =  WHJ_data[3];
//       T WHJ01       =  WHJ_data[1];
//       T WHJ10       =  WHJ_data[2];

//       out_data[0]   =  H00*WHJ11;
//       out_data[1]   =  H01*WHJ11;
//       out_data[2]   =  H10*WHJ11;
//       out_data[3]   =  H11*WHJ11;
//       out_data[4]   =  -H00*WHJ10;
//       out_data[5]   =  -H01*WHJ10;
//       out_data[6]   =  -H10*WHJ10;
//       out_data[7]   =  -H11*WHJ10;
//       out_data[8]   =  -H00*WHJ01;
//       out_data[9]   =  -H01*WHJ01;
//       out_data[10]  =  -H10*WHJ01;
//       out_data[11]  =  -H11*WHJ01;
//       out_data[12]  =  H00*WHJ00;
//       out_data[13]  =  H01*WHJ00;
//       out_data[14]  =  H10*WHJ00;
//       out_data[15]  =  H11*WHJ00;
//    }
//}

//tensor<real,3> H_HJ(const tensor<real,3> &F,const tensor<real,3> &H, const tensor<real,3> &WHJ,  const int &ndim) {
//    //--------------------------------------------------------------
//    // Obtain FH contribution in the Hessian for every Gauss point
//    //--------------------------------------------------------------
//    auto ngauss        =  F.dimension(2);
//    tensor<real,5> out(ndim,ndim,ndim,ndim,ngauss);
//    for (auto i=0; i<ngauss; ++i) {
//        tensor<real,2> F_2d   =  F.chip(i,2);
//        tensor<real,2> H_2d   =  H.chip(i,2);
//        tensor<real,2> WHJ_2d =  WHJ.chip(i,2);
//        __Hessian_HJ_contribution__(F_2d.data(),H_2d.data(),WHJ_2d.data(), out.data()+ndim*ndim*ndim*ndim*i,ndim);
//    }
//    return out;
//}

//template<typename T>
//SC_INLINE void __Hessian_JF_contribution__(const T *__restrict__ a_data, const T *__restrict__ b_data, T *__restrict__ out_data, const int &ndim) {
//    //--------------------------------------------------------------
//    // Contribution JF in the Hessian operator
//    //--------------------------------------------------------------
//    if (ndim==3){
//       T H00         =  a_data[0];
//       T H11         =  a_data[4];
//       T H22         =  a_data[8];
//       T H01         =  a_data[1];
//       T H02         =  a_data[2];
//       T H12         =  a_data[5];
//       T H10         =  a_data[3];
//       T H20         =  a_data[6];
//       T H21         =  a_data[7];

//       T WJF00       =  b_data[0];
//       T WJF11       =  b_data[4];
//       T WJF22       =  b_data[8];
//       T WJF01       =  b_data[1];
//       T WJF02       =  b_data[2];
//       T WJF12       =  b_data[5];
//       T WJF10       =  b_data[3];
//       T WJF20       =  b_data[6];
//       T WJF21       =  b_data[7];

//       out_data[0]   =  H00*WJF00;
//       out_data[1]   =  H01*WJF00;
//       out_data[2]   =  H02*WJF00;
//       out_data[3]   =  H10*WJF00;
//       out_data[4]   =  H11*WJF00;
//       out_data[5]   =  H12*WJF00;
//       out_data[6]   =  H20*WJF00;
//       out_data[7]   =  H21*WJF00;
//       out_data[8]   =  H22*WJF00;
//       out_data[9]   =  H00*WJF01;
//       out_data[10]  =  H01*WJF01;
//       out_data[11]  =  H02*WJF01;
//       out_data[12]  =  H10*WJF01;
//       out_data[13]  =  H11*WJF01;
//       out_data[14]  =  H12*WJF01;
//       out_data[15]  =  H20*WJF01;
//       out_data[16]  =  H21*WJF01;
//       out_data[17]  =  H22*WJF01;
//       out_data[18]  =  H00*WJF02;
//       out_data[19]  =  H01*WJF02;
//       out_data[20]  =  H02*WJF02;
//       out_data[21]  =  H10*WJF02;
//       out_data[22]  =  H11*WJF02;
//       out_data[23]  =  H12*WJF02;
//       out_data[24]  =  H20*WJF02;
//       out_data[25]  =  H21*WJF02;
//       out_data[26]  =  H22*WJF02;
//       out_data[27]  =  H00*WJF10;
//       out_data[28]  =  H01*WJF10;
//       out_data[29]  =  H02*WJF10;
//       out_data[30]  =  H10*WJF10;
//       out_data[31]  =  H11*WJF10;
//       out_data[32]  =  H12*WJF10;
//       out_data[33]  =  H20*WJF10;
//       out_data[34]  =  H21*WJF10;
//       out_data[35]  =  H22*WJF10;
//       out_data[36]  =  H00*WJF11;
//       out_data[37]  =  H01*WJF11;
//       out_data[38]  =  H02*WJF11;
//       out_data[39]  =  H10*WJF11;
//       out_data[40]  =  H11*WJF11;
//       out_data[41]  =  H12*WJF11;
//       out_data[42]  =  H20*WJF11;
//       out_data[43]  =  H21*WJF11;
//       out_data[44]  =  H22*WJF11;
//       out_data[45]  =  H00*WJF12;
//       out_data[46]  =  H01*WJF12;
//       out_data[47]  =  H02*WJF12;
//       out_data[48]  =  H10*WJF12;
//       out_data[49]  =  H11*WJF12;
//       out_data[50]  =  H12*WJF12;
//       out_data[51]  =  H20*WJF12;
//       out_data[52]  =  H21*WJF12;
//       out_data[53]  =  H22*WJF12;
//       out_data[54]  =  H00*WJF20;
//       out_data[55]  =  H01*WJF20;
//       out_data[56]  =  H02*WJF20;
//       out_data[57]  =  H10*WJF20;
//       out_data[58]  =  H11*WJF20;
//       out_data[59]  =  H12*WJF20;
//       out_data[60]  =  H20*WJF20;
//       out_data[61]  =  H21*WJF20;
//       out_data[62]  =  H22*WJF20;
//       out_data[63]  =  H00*WJF21;
//       out_data[64]  =  H01*WJF21;
//       out_data[65]  =  H02*WJF21;
//       out_data[66]  =  H10*WJF21;
//       out_data[67]  =  H11*WJF21;
//       out_data[68]  =  H12*WJF21;
//       out_data[69]  =  H20*WJF21;
//       out_data[70]  =  H21*WJF21;
//       out_data[71]  =  H22*WJF21;
//       out_data[72]  =  H00*WJF22;
//       out_data[73]  =  H01*WJF22;
//       out_data[74]  =  H02*WJF22;
//       out_data[75]  =  H10*WJF22;
//       out_data[76]  =  H11*WJF22;
//       out_data[77]  =  H12*WJF22;
//       out_data[78]  =  H20*WJF22;
//       out_data[79]  =  H21*WJF22;
//       out_data[80]  =  H22*WJF22;
//    }
//    else{
//        T H00         =  a_data[0];
//        T H11         =  a_data[3];
//        T H01         =  a_data[1];
//        T H10         =  a_data[2];

//        T WJF00       =  b_data[0];
//        T WJF11       =  b_data[3];
//        T WJF01       =  b_data[1];
//        T WJF10       =  b_data[2];

//        out_data[0]   =  H00*WJF00;
//        out_data[1]   =  H01*WJF00;
//        out_data[2]   =  H10*WJF00;
//        out_data[3]   =  H11*WJF00;
//        out_data[4]   =  H00*WJF01;
//        out_data[5]   =  H01*WJF01;
//        out_data[6]   =  H10*WJF01;
//        out_data[7]   =  H11*WJF01;
//        out_data[8]   =  H00*WJF10;
//        out_data[9]   =  H01*WJF10;
//        out_data[10]  =  H10*WJF10;
//        out_data[11]  =  H11*WJF10;
//        out_data[12]  =  H00*WJF11;
//        out_data[13]  =  H01*WJF11;
//        out_data[14]  =  H10*WJF11;
//        out_data[15]  =  H11*WJF11;     }
//}

//tensor<real,5> H_JF(const tensor<real,3> &F,const tensor<real,3> &WJF,  const int &ndim) {
//    //--------------------------------------------------------------
//    // Obtain FH contribution in the Hessian for every Gauss point
//    //--------------------------------------------------------------
//    auto ngauss        =  F.dimension(2);
//    tensor<real,5> out(ndim,ndim,ndim,ndim,ngauss);
//    for (auto i=0; i<ngauss; ++i) {
//        tensor<real,2> F_2d   =  F.chip(i,2);
//        tensor<real,2> WJF_2d =  WJF.chip(i,2);
//        __Hessian_JF_contribution__(F_2d.data(),WJF_2d.data(), out.data()+ndim*ndim*ndim*ndim*i,ndim);
//    }
//    return out;
//}

//template<typename T>
//SC_INLINE void __Hessian_JH_contribution__(const T *__restrict__ F_data, const T *__restrict__ H_data, const T *__restrict__ WJH_data, T *__restrict__ out_data, const int &ndim) {
//    //--------------------------------------------------------------
//    // Contribution HF in the Hessian operator
//    //--------------------------------------------------------------
//    if (ndim==3){
//       T F00         =  F_data[0];
//       T F11         =  F_data[4];
//       T F22         =  F_data[8];
//       T F01         =  F_data[1];
//       T F02         =  F_data[2];
//       T F12         =  F_data[5];
//       T F10         =  F_data[3];
//       T F20         =  F_data[6];
//       T F21         =  F_data[7];

//       T H00         =  H_data[0];
//       T H11         =  H_data[4];
//       T H22         =  H_data[8];
//       T H01         =  H_data[1];
//       T H02         =  H_data[2];
//       T H12         =  H_data[5];
//       T H10         =  H_data[3];
//       T H20         =  H_data[6];
//       T H21         =  H_data[7];

//       T WJH00       =  WJH_data[0];
//       T WJH11       =  WJH_data[4];
//       T WJH22       =  WJH_data[8];
//       T WJH01       =  WJH_data[1];
//       T WJH02       =  WJH_data[2];
//       T WJH12       =  WJH_data[5];
//       T WJH10       =  WJH_data[3];
//       T WJH20       =  WJH_data[6];
//       T WJH21       =  WJH_data[7];

//       out_data[0]   =  F11*H00*WJH22 - F12*H00*WJH21 - F21*H00*WJH12 + F22*H00*WJH11;
//       out_data[1]   =  F12*H00*WJH20 - F10*H00*WJH22 + F20*H00*WJH12 - F22*H00*WJH10;
//       out_data[2]   =  F10*H00*WJH21 - F11*H00*WJH20 - F20*H00*WJH11 + F21*H00*WJH10;
//       out_data[3]   =  F02*H00*WJH21 - F01*H00*WJH22 + F21*H00*WJH02 - F22*H00*WJH01;
//       out_data[4]   =  F00*H00*WJH22 - F02*H00*WJH20 - F20*H00*WJH02 + F22*H00*WJH00;
//       out_data[5]   =  F01*H00*WJH20 - F00*H00*WJH21 + F20*H00*WJH01 - F21*H00*WJH00;
//       out_data[6]   =  F01*H00*WJH12 - F02*H00*WJH11 - F11*H00*WJH02 + F12*H00*WJH01;
//       out_data[7]   =  F02*H00*WJH10 - F00*H00*WJH12 + F10*H00*WJH02 - F12*H00*WJH00;
//       out_data[8]   =  F00*H00*WJH11 - F01*H00*WJH10 - F10*H00*WJH01 + F11*H00*WJH00;
//       out_data[9]   =  F11*H01*WJH22 - F12*H01*WJH21 - F21*H01*WJH12 + F22*H01*WJH11;
//       out_data[10]  =  F12*H01*WJH20 - F10*H01*WJH22 + F20*H01*WJH12 - F22*H01*WJH10;
//       out_data[11]  =  F10*H01*WJH21 - F11*H01*WJH20 - F20*H01*WJH11 + F21*H01*WJH10;
//       out_data[12]  =  F02*H01*WJH21 - F01*H01*WJH22 + F21*H01*WJH02 - F22*H01*WJH01;
//       out_data[13]  =  F00*H01*WJH22 - F02*H01*WJH20 - F20*H01*WJH02 + F22*H01*WJH00;
//       out_data[14]  =  F01*H01*WJH20 - F00*H01*WJH21 + F20*H01*WJH01 - F21*H01*WJH00;
//       out_data[15]  =  F01*H01*WJH12 - F02*H01*WJH11 - F11*H01*WJH02 + F12*H01*WJH01;
//       out_data[16]  =  F02*H01*WJH10 - F00*H01*WJH12 + F10*H01*WJH02 - F12*H01*WJH00;
//       out_data[17]  =  F00*H01*WJH11 - F01*H01*WJH10 - F10*H01*WJH01 + F11*H01*WJH00;
//       out_data[18]  =  F11*H02*WJH22 - F12*H02*WJH21 - F21*H02*WJH12 + F22*H02*WJH11;
//       out_data[19]  =  F12*H02*WJH20 - F10*H02*WJH22 + F20*H02*WJH12 - F22*H02*WJH10;
//       out_data[20]  =  F10*H02*WJH21 - F11*H02*WJH20 - F20*H02*WJH11 + F21*H02*WJH10;
//       out_data[21]  =  F02*H02*WJH21 - F01*H02*WJH22 + F21*H02*WJH02 - F22*H02*WJH01;
//       out_data[22]  =  F00*H02*WJH22 - F02*H02*WJH20 - F20*H02*WJH02 + F22*H02*WJH00;
//       out_data[23]  =  F01*H02*WJH20 - F00*H02*WJH21 + F20*H02*WJH01 - F21*H02*WJH00;
//       out_data[24]  =  F01*H02*WJH12 - F02*H02*WJH11 - F11*H02*WJH02 + F12*H02*WJH01;
//       out_data[25]  =  F02*H02*WJH10 - F00*H02*WJH12 + F10*H02*WJH02 - F12*H02*WJH00;
//       out_data[26]  =  F00*H02*WJH11 - F01*H02*WJH10 - F10*H02*WJH01 + F11*H02*WJH00;
//       out_data[27]  =  F11*H10*WJH22 - F12*H10*WJH21 - F21*H10*WJH12 + F22*H10*WJH11;
//       out_data[28]  =  F12*H10*WJH20 - F10*H10*WJH22 + F20*H10*WJH12 - F22*H10*WJH10;
//       out_data[29]  =  F10*H10*WJH21 - F11*H10*WJH20 - F20*H10*WJH11 + F21*H10*WJH10;
//       out_data[30]  =  F02*H10*WJH21 - F01*H10*WJH22 + F21*H10*WJH02 - F22*H10*WJH01;
//       out_data[31]  =  F00*H10*WJH22 - F02*H10*WJH20 - F20*H10*WJH02 + F22*H10*WJH00;
//       out_data[32]  =  F01*H10*WJH20 - F00*H10*WJH21 + F20*H10*WJH01 - F21*H10*WJH00;
//       out_data[33]  =  F01*H10*WJH12 - F02*H10*WJH11 - F11*H10*WJH02 + F12*H10*WJH01;
//       out_data[34]  =  F02*H10*WJH10 - F00*H10*WJH12 + F10*H10*WJH02 - F12*H10*WJH00;
//       out_data[35]  =  F00*H10*WJH11 - F01*H10*WJH10 - F10*H10*WJH01 + F11*H10*WJH00;
//       out_data[36]  =  F11*H11*WJH22 - F12*H11*WJH21 - F21*H11*WJH12 + F22*H11*WJH11;
//       out_data[37]  =  F12*H11*WJH20 - F10*H11*WJH22 + F20*H11*WJH12 - F22*H11*WJH10;
//       out_data[38]  =  F10*H11*WJH21 - F11*H11*WJH20 - F20*H11*WJH11 + F21*H11*WJH10;
//       out_data[39]  =  F02*H11*WJH21 - F01*H11*WJH22 + F21*H11*WJH02 - F22*H11*WJH01;
//       out_data[40]  =  F00*H11*WJH22 - F02*H11*WJH20 - F20*H11*WJH02 + F22*H11*WJH00;
//       out_data[41]  =  F01*H11*WJH20 - F00*H11*WJH21 + F20*H11*WJH01 - F21*H11*WJH00;
//       out_data[42]  =  F01*H11*WJH12 - F02*H11*WJH11 - F11*H11*WJH02 + F12*H11*WJH01;
//       out_data[43]  =  F02*H11*WJH10 - F00*H11*WJH12 + F10*H11*WJH02 - F12*H11*WJH00;
//       out_data[44]  =  F00*H11*WJH11 - F01*H11*WJH10 - F10*H11*WJH01 + F11*H11*WJH00;
//       out_data[45]  =  F11*H12*WJH22 - F12*H12*WJH21 - F21*H12*WJH12 + F22*H12*WJH11;
//       out_data[46]  =  F12*H12*WJH20 - F10*H12*WJH22 + F20*H12*WJH12 - F22*H12*WJH10;
//       out_data[47]  =  F10*H12*WJH21 - F11*H12*WJH20 - F20*H12*WJH11 + F21*H12*WJH10;
//       out_data[48]  =  F02*H12*WJH21 - F01*H12*WJH22 + F21*H12*WJH02 - F22*H12*WJH01;
//       out_data[49]  =  F00*H12*WJH22 - F02*H12*WJH20 - F20*H12*WJH02 + F22*H12*WJH00;
//       out_data[50]  =  F01*H12*WJH20 - F00*H12*WJH21 + F20*H12*WJH01 - F21*H12*WJH00;
//       out_data[51]  =  F01*H12*WJH12 - F02*H12*WJH11 - F11*H12*WJH02 + F12*H12*WJH01;
//       out_data[52]  =  F02*H12*WJH10 - F00*H12*WJH12 + F10*H12*WJH02 - F12*H12*WJH00;
//       out_data[53]  =  F00*H12*WJH11 - F01*H12*WJH10 - F10*H12*WJH01 + F11*H12*WJH00;
//       out_data[54]  =  F11*H20*WJH22 - F12*H20*WJH21 - F21*H20*WJH12 + F22*H20*WJH11;
//       out_data[55]  =  F12*H20*WJH20 - F10*H20*WJH22 + F20*H20*WJH12 - F22*H20*WJH10;
//       out_data[56]  =  F10*H20*WJH21 - F11*H20*WJH20 - F20*H20*WJH11 + F21*H20*WJH10;
//       out_data[57]  =  F02*H20*WJH21 - F01*H20*WJH22 + F21*H20*WJH02 - F22*H20*WJH01;
//       out_data[58]  =  F00*H20*WJH22 - F02*H20*WJH20 - F20*H20*WJH02 + F22*H20*WJH00;
//       out_data[59]  =  F01*H20*WJH20 - F00*H20*WJH21 + F20*H20*WJH01 - F21*H20*WJH00;
//       out_data[60]  =  F01*H20*WJH12 - F02*H20*WJH11 - F11*H20*WJH02 + F12*H20*WJH01;
//       out_data[61]  =  F02*H20*WJH10 - F00*H20*WJH12 + F10*H20*WJH02 - F12*H20*WJH00;
//       out_data[62]  =  F00*H20*WJH11 - F01*H20*WJH10 - F10*H20*WJH01 + F11*H20*WJH00;
//       out_data[63]  =  F11*H21*WJH22 - F12*H21*WJH21 - F21*H21*WJH12 + F22*H21*WJH11;
//       out_data[64]  =  F12*H21*WJH20 - F10*H21*WJH22 + F20*H21*WJH12 - F22*H21*WJH10;
//       out_data[65]  =  F10*H21*WJH21 - F11*H21*WJH20 - F20*H21*WJH11 + F21*H21*WJH10;
//       out_data[66]  =  F02*H21*WJH21 - F01*H21*WJH22 + F21*H21*WJH02 - F22*H21*WJH01;
//       out_data[67]  =  F00*H21*WJH22 - F02*H21*WJH20 - F20*H21*WJH02 + F22*H21*WJH00;
//       out_data[68]  =  F01*H21*WJH20 - F00*H21*WJH21 + F20*H21*WJH01 - F21*H21*WJH00;
//       out_data[69]  =  F01*H21*WJH12 - F02*H21*WJH11 - F11*H21*WJH02 + F12*H21*WJH01;
//       out_data[70]  =  F02*H21*WJH10 - F00*H21*WJH12 + F10*H21*WJH02 - F12*H21*WJH00;
//       out_data[71]  =  F00*H21*WJH11 - F01*H21*WJH10 - F10*H21*WJH01 + F11*H21*WJH00;
//       out_data[72]  =  F11*H22*WJH22 - F12*H22*WJH21 - F21*H22*WJH12 + F22*H22*WJH11;
//       out_data[73]  =  F12*H22*WJH20 - F10*H22*WJH22 + F20*H22*WJH12 - F22*H22*WJH10;
//       out_data[74]  =  F10*H22*WJH21 - F11*H22*WJH20 - F20*H22*WJH11 + F21*H22*WJH10;
//       out_data[75]  =  F02*H22*WJH21 - F01*H22*WJH22 + F21*H22*WJH02 - F22*H22*WJH01;
//       out_data[76]  =  F00*H22*WJH22 - F02*H22*WJH20 - F20*H22*WJH02 + F22*H22*WJH00;
//       out_data[77]  =  F01*H22*WJH20 - F00*H22*WJH21 + F20*H22*WJH01 - F21*H22*WJH00;
//       out_data[78]  =  F01*H22*WJH12 - F02*H22*WJH11 - F11*H22*WJH02 + F12*H22*WJH01;
//       out_data[79]  =  F02*H22*WJH10 - F00*H22*WJH12 + F10*H22*WJH02 - F12*H22*WJH00;
//       out_data[80]  =  F00*H22*WJH11 - F01*H22*WJH10 - F10*H22*WJH01 + F11*H22*WJH00;     }
//    else{
//       T H00         =  H_data[0];
//       T H11         =  H_data[3];
//       T H01         =  H_data[1];
//       T H10         =  H_data[2];

//       T WJH00       =  WJH_data[0];
//       T WJH11       =  WJH_data[3];
//       T WJH01       =  WJH_data[1];
//       T WJH10       =  WJH_data[2];

//       out_data[0]   =  H00*WJH11;
//       out_data[1]   =  -H00*WJH10;
//       out_data[2]   =  -H00*WJH01;
//       out_data[3]   =  H00*WJH00;
//       out_data[4]   =  H01*WJH11;
//       out_data[5]   =  -H01*WJH10;
//       out_data[6]   =  -H01*WJH01;
//       out_data[7]   =  H01*WJH00;
//       out_data[8]   =  H10*WJH11;
//       out_data[9]   =  -H10*WJH10;
//       out_data[10]  =  -H10*WJH01;
//       out_data[11]  =  H10*WJH00;
//       out_data[12]  =  H11*WJH11;
//       out_data[13]  =  -H11*WJH10;
//       out_data[14]  =  -H11*WJH01;
//       out_data[15]  =  H11*WJH00;     }
//}

//tensor<real,5> H_JH(const tensor<real,3> &F,const tensor<real,3> &H, const tensor<real,3> &WJH,  const int &ndim) {
//    //--------------------------------------------------------------
//    // Obtain FH contribution in the Hessian for every Gauss point
//    //--------------------------------------------------------------
//    auto ngauss        =  F.dimension(2);
//    tensor<real,5> out(ndim,ndim,ndim,ndim,ngauss);
//    for (auto i=0; i<ngauss; ++i) {
//        tensor<real,2> F_2d   =  F.chip(i,2);
//        tensor<real,2> H_2d   =  H.chip(i,2);
//        tensor<real,2> WJH_2d =  WJH.chip(i,2);
//        __Hessian_JH_contribution__(F_2d.data(),H_2d.data(),WJH_2d.data(), out.data()+ndim*ndim*ndim*ndim*i,ndim);
//    }
//    return out;
//}

template<typename T>
SC_INLINE void __Hessian_JJ_contribution__(const T *__restrict__ a_data, const real WJJ, T *__restrict__ out_data, const int &ndim) {
    //--------------------------------------------------------------
    // Contribution JJ in the Hessian operator
    //--------------------------------------------------------------
    if (ndim==3){
       T H00         =  a_data[0];
       T H11         =  a_data[4];
       T H22         =  a_data[8];
       T H01         =  a_data[1];
       T H02         =  a_data[2];
       T H12         =  a_data[5];
       T H10         =  a_data[3];
       T H20         =  a_data[6];
       T H21         =  a_data[7];

       out_data[0]   =  H00*H00*WJJ;
       out_data[1]   =  H00*H01*WJJ;
       out_data[2]   =  H00*H02*WJJ;
       out_data[3]   =  H00*H10*WJJ;
       out_data[4]   =  H00*H11*WJJ;
       out_data[5]   =  H00*H12*WJJ;
       out_data[6]   =  H00*H20*WJJ;
       out_data[7]   =  H00*H21*WJJ;
       out_data[8]   =  H00*H22*WJJ;
       out_data[9]   =  H00*H01*WJJ;
       out_data[10]  =  H01*H01*WJJ;
       out_data[11]  =  H01*H02*WJJ;
       out_data[12]  =  H01*H10*WJJ;
       out_data[13]  =  H01*H11*WJJ;
       out_data[14]  =  H01*H12*WJJ;
       out_data[15]  =  H01*H20*WJJ;
       out_data[16]  =  H01*H21*WJJ;
       out_data[17]  =  H01*H22*WJJ;
       out_data[18]  =  H00*H02*WJJ;
       out_data[19]  =  H01*H02*WJJ;
       out_data[20]  =  H02*H02*WJJ;
       out_data[21]  =  H02*H10*WJJ;
       out_data[22]  =  H02*H11*WJJ;
       out_data[23]  =  H02*H12*WJJ;
       out_data[24]  =  H02*H20*WJJ;
       out_data[25]  =  H02*H21*WJJ;
       out_data[26]  =  H02*H22*WJJ;
       out_data[27]  =  H00*H10*WJJ;
       out_data[28]  =  H01*H10*WJJ;
       out_data[29]  =  H02*H10*WJJ;
       out_data[30]  =  H10*H10*WJJ;
       out_data[31]  =  H10*H11*WJJ;
       out_data[32]  =  H10*H12*WJJ;
       out_data[33]  =  H10*H20*WJJ;
       out_data[34]  =  H10*H21*WJJ;
       out_data[35]  =  H10*H22*WJJ;
       out_data[36]  =  H00*H11*WJJ;
       out_data[37]  =  H01*H11*WJJ;
       out_data[38]  =  H02*H11*WJJ;
       out_data[39]  =  H10*H11*WJJ;
       out_data[40]  =  H11*H11*WJJ;
       out_data[41]  =  H11*H12*WJJ;
       out_data[42]  =  H11*H20*WJJ;
       out_data[43]  =  H11*H21*WJJ;
       out_data[44]  =  H11*H22*WJJ;
       out_data[45]  =  H00*H12*WJJ;
       out_data[46]  =  H01*H12*WJJ;
       out_data[47]  =  H02*H12*WJJ;
       out_data[48]  =  H10*H12*WJJ;
       out_data[49]  =  H11*H12*WJJ;
       out_data[50]  =  H12*H12*WJJ;
       out_data[51]  =  H12*H20*WJJ;
       out_data[52]  =  H12*H21*WJJ;
       out_data[53]  =  H12*H22*WJJ;
       out_data[54]  =  H00*H20*WJJ;
       out_data[55]  =  H01*H20*WJJ;
       out_data[56]  =  H02*H20*WJJ;
       out_data[57]  =  H10*H20*WJJ;
       out_data[58]  =  H11*H20*WJJ;
       out_data[59]  =  H12*H20*WJJ;
       out_data[60]  =  H20*H20*WJJ;
       out_data[61]  =  H20*H21*WJJ;
       out_data[62]  =  H20*H22*WJJ;
       out_data[63]  =  H00*H21*WJJ;
       out_data[64]  =  H01*H21*WJJ;
       out_data[65]  =  H02*H21*WJJ;
       out_data[66]  =  H10*H21*WJJ;
       out_data[67]  =  H11*H21*WJJ;
       out_data[68]  =  H12*H21*WJJ;
       out_data[69]  =  H20*H21*WJJ;
       out_data[70]  =  H21*H21*WJJ;
       out_data[71]  =  H21*H22*WJJ;
       out_data[72]  =  H00*H22*WJJ;
       out_data[73]  =  H01*H22*WJJ;
       out_data[74]  =  H02*H22*WJJ;
       out_data[75]  =  H10*H22*WJJ;
       out_data[76]  =  H11*H22*WJJ;
       out_data[77]  =  H12*H22*WJJ;
       out_data[78]  =  H20*H22*WJJ;
       out_data[79]  =  H21*H22*WJJ;
       out_data[80]  =  H22*H22*WJJ;
    }
    else{
       T H00         =  a_data[0];
       T H01         =  a_data[1];
       T H10         =  a_data[2];
       T H11         =  a_data[3];

       out_data[0]   =  H00*H00*WJJ;
       out_data[1]   =  H00*H01*WJJ;
       out_data[2]   =  H00*H10*WJJ;
       out_data[3]   =  H00*H11*WJJ;
       out_data[4]   =  H00*H01*WJJ;
       out_data[5]   =  H01*H01*WJJ;
       out_data[6]   =  H01*H10*WJJ;
       out_data[7]   =  H01*H11*WJJ;
       out_data[8]   =  H00*H10*WJJ;
       out_data[9]   =  H01*H10*WJJ;
       out_data[10]  =  H10*H10*WJJ;
       out_data[11]  =  H10*H11*WJJ;
       out_data[12]  =  H00*H11*WJJ;
       out_data[13]  =  H01*H11*WJJ;
       out_data[14]  =  H10*H11*WJJ;
       out_data[15]  =  H11*H11*WJJ;
    }
}

tensor<real,5> H_JJ(const tensor<real,3> &H,const tensor<real,1> &WJJ,  const int &ndim) {
    //--------------------------------------------------------------
    // Obtain JJ contribution in the Hessian for every Gauss point
    //--------------------------------------------------------------
    auto ngauss        =  H.dimension(0);
    tensor<real,5> out(ngauss,ndim,ndim,ndim,ndim);
    for (auto i=0; i<ngauss; ++i) {
        tensor<real,2> H_2d   =  H.chip(i,0);
        real WJJ_2d =  WJJ(i);
        __Hessian_JJ_contribution__(H_2d.data(),WJJ_2d, out.data()+ndim*ndim*ndim*ndim*i,ndim);
    }
    return out;
}

template<typename T>
SC_INLINE void __Hessian_geom_contribution__(const T *__restrict__ F_data, const T *__restrict__ sigmaH_data, const real SigmaJ, T *__restrict__ out_data, const int &ndim) {
    //--------------------------------------------------------------
    // Geometrical contribution in the Hessian operator
    //--------------------------------------------------------------
    if (ndim==3) {
       T F00         =  F_data[0];
       T F01         =  F_data[1];
       T F02         =  F_data[2];
       T F10         =  F_data[3];
       T F11         =  F_data[4];
       T F12         =  F_data[5];
       T F20         =  F_data[6];
       T F21         =  F_data[7];
       T F22         =  F_data[8];

       T sigma_H00   =  sigmaH_data[0];
       T sigma_H01   =  sigmaH_data[1];
       T sigma_H02   =  sigmaH_data[2];
       T sigma_H10   =  sigmaH_data[3];
       T sigma_H11   =  sigmaH_data[4];
       T sigma_H12   =  sigmaH_data[5];
       T sigma_H20   =  sigmaH_data[6];
       T sigma_H21   =  sigmaH_data[7];
       T sigma_H22   =  sigmaH_data[8];

       out_data[0]   =  0.0;
       out_data[1]   =  0.0;
       out_data[2]   =  0.0;
       out_data[3]   =  0.0;
       out_data[4]   =  sigma_H22 + F22*SigmaJ;
       out_data[5]   =  - sigma_H21 - F21*SigmaJ;
       out_data[6]   =  0.0;
       out_data[7]   =  - sigma_H12 - F12*SigmaJ;
       out_data[8]   =  sigma_H11 + F11*SigmaJ;
       out_data[9]   =  0.0;
       out_data[10]  =  0.0;
       out_data[11]  =  0.0;
       out_data[12]  =  - sigma_H22 - F22*SigmaJ;
       out_data[13]  =  0.0;
       out_data[14]  =  sigma_H20 + F20*SigmaJ;
       out_data[15]  =  sigma_H12 + F12*SigmaJ;
       out_data[16]  =  0.0;
       out_data[17]  =  - sigma_H10 - F10*SigmaJ;
       out_data[18]  =  0.0;
       out_data[19]  =  0.0;
       out_data[20]  =  0.0;
       out_data[21]  =  sigma_H21 + F21*SigmaJ;
       out_data[22]  =  - sigma_H20 - F20*SigmaJ;
       out_data[23]  =  0.0;
       out_data[24]  =  - sigma_H11 - F11*SigmaJ;
       out_data[25]  =  sigma_H10 + F10*SigmaJ;
       out_data[26]  =  0.0;
       out_data[27]  =  0.0;
       out_data[28]  =  - sigma_H22 - F22*SigmaJ;
       out_data[29]  =  sigma_H21 + F21*SigmaJ;
       out_data[30]  =  0.0;
       out_data[31]  =  0.0;
       out_data[32]  =  0.0;
       out_data[33]  =  0.0;
       out_data[34]  =  sigma_H02 + F02*SigmaJ;
       out_data[35]  =  - sigma_H01 - F01*SigmaJ;
       out_data[36]  =  sigma_H22 + F22*SigmaJ;
       out_data[37]  =  0.0;
       out_data[38]  =  - sigma_H20 - F20*SigmaJ;
       out_data[39]  =  0.0;
       out_data[40]  =  0.0;
       out_data[41]  =  0.0;
       out_data[42]  =  - sigma_H02 - F02*SigmaJ;
       out_data[43]  =  0.0;
       out_data[44]  =  sigma_H00 + F00*SigmaJ;
       out_data[45]  =  - sigma_H21 - F21*SigmaJ;
       out_data[46]  =  sigma_H20 + F20*SigmaJ;
       out_data[47]  =  0.0;
       out_data[48]  =  0.0;
       out_data[49]  =  0.0;
       out_data[50]  =  0.0;
       out_data[51]  =  sigma_H01 + F01*SigmaJ;
       out_data[52]  =  - sigma_H00 - F00*SigmaJ;
       out_data[53]  =  0.0;
       out_data[54]  =  0.0;
       out_data[55]  =  sigma_H12 + F12*SigmaJ;
       out_data[56]  =  - sigma_H11 - F11*SigmaJ;
       out_data[57]  =  0.0;
       out_data[58]  =  - sigma_H02 - F02*SigmaJ;
       out_data[59]  =  sigma_H01 + F01*SigmaJ;
       out_data[60]  =  0.0;
       out_data[61]  =  0.0;
       out_data[62]  =  0.0;
       out_data[63]  =  - sigma_H12 - F12*SigmaJ;
       out_data[64]  =  0.0;
       out_data[65]  =  sigma_H10 + F10*SigmaJ;
       out_data[66]  =  sigma_H02 + F02*SigmaJ;
       out_data[67]  =  0.0;
       out_data[68]  =  - sigma_H00 - F00*SigmaJ;
       out_data[69]  =  0.0;
       out_data[70]  =  0.0;
       out_data[71]  =  0.0;
       out_data[72]  =  sigma_H11 + F11*SigmaJ;
       out_data[73]  =  - sigma_H10 - F10*SigmaJ;
       out_data[74]  =  0.0;
       out_data[75]  =  - sigma_H01 - F01*SigmaJ;
       out_data[76]  =  sigma_H00 + F00*SigmaJ;
       out_data[77]  =  0.0;
       out_data[78]  =  0.0;
       out_data[79]  =  0.0;
       out_data[80]  =  0.0;
       }
   else {
       out_data[0]   =  0.0;
       out_data[1]   =  0.0;
       out_data[2]   =  0.0;
       out_data[3]   =  SigmaJ;
       out_data[4]   =  0.0;
       out_data[5]   =  0.0;
       out_data[6]   =  -SigmaJ;
       out_data[7]   =  0.0;
       out_data[8]   =  0.0;
       out_data[9]   =  -SigmaJ;
       out_data[10]  =  0.0;
       out_data[11]  =  0.0;
       out_data[12]  =  SigmaJ;
       out_data[13]  =  0.0;
       out_data[14]  =  0.0;
       out_data[15]  =  0.0;
   }
}

tensor<real,5> H_geom(const tensor<real,3> &F, const tensor<real,3> &sigma_H, const tensor<real,1> &sigma_J,  const int &ndim) {
    //--------------------------------------------------------------
    // Obtain JJ contribution in the Hessian for every Gauss point
    //--------------------------------------------------------------
    auto ngauss        =  F.dimension(0);
    tensor<real,5> out(ngauss,ndim,ndim,ndim,ndim);
    for (auto i=0; i<ngauss; ++i) {
        tensor<real,2> F_2d         =  F.chip(i,0);
        tensor<real,2> sigma_H_2d   =  sigma_H.chip(i,0);
        real sigma_J_2d             =  sigma_J(i);
        __Hessian_geom_contribution__(F_2d.data(),sigma_H_2d.data(),sigma_J_2d, out.data()+ndim*ndim*ndim*ndim*i,ndim);
    }
    return out;
}






}
#endif // BACKEND_CONSTITUTIVETENSORS_H
