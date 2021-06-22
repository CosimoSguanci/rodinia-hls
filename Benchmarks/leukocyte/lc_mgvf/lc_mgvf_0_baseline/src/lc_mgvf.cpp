#include"lc_mgvf.h"
#include <ap_fixed.h>
#include <hls_math.h>
#include <cmath>
#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>


extern "C" {

ap_fixed<8, 1> heaviside(ap_fixed<8, 1> x) {
    // A simpler, faster approximation of the Heaviside function
    ap_fixed<8, 1> out = 0.0;
    if (x > -0.0001) out = 0.5;
    if (x >  0.0001) out = 1.0;
    return out; 
}

float lc_mgvf(ap_fixed<8, 1> result[GRID_ROWS * GRID_COLS], ap_fixed<8, 1> imgvf[GRID_ROWS * GRID_COLS], ap_fixed<8, 1> I[GRID_ROWS * GRID_COLS])
{

    ap_fixed<8, 1> total_diff = 0.0;


    for (int i = 0; i < GRID_ROWS; i++) { 
        for (int j = 0; j < GRID_COLS; j++) {

            ap_fixed<8, 1> old_val = imgvf[i * GRID_COLS + j];

            ap_fixed<8, 1> UL;
            if(i == 0              ||  j == 0              ) {
                UL = 0;
            }
            else {
                UL = imgvf[(i - 1   ) * GRID_COLS + (j - 1  )] - old_val;
            }

            
            ap_fixed<8, 1> U;
            if(i == 0                                      ) {
                U = 0;
            }
            else {
                U = imgvf[(i - 1   ) * GRID_COLS + (j      )] - old_val;
            }

            ap_fixed<8, 1> UR;
            if(i == 0              ||  j == GRID_COLS - 1  ) {
                UR = 0;
            }
            else {
                UR = imgvf[(i - 1   ) * GRID_COLS + (j + 1  )] - old_val;
            }    
            

            ap_fixed<8, 1> L;
            if(                        j == 0              ) {
                L = 0;
            }
            else {
                L = imgvf[(i       ) * GRID_COLS + (j - 1  )] - old_val;
            }

            
            ap_fixed<8, 1> R;
            if(                        j == GRID_COLS - 1  ) {
                R = 0;
            }
            else {
                R = imgvf[(i       ) * GRID_COLS + (j + 1  )] - old_val;
            }

            ap_fixed<8, 1> DL;
            if(i == GRID_ROWS - 1  ||  j == 0              ) {
                DL = 0;
            }
            else {
                DL = imgvf[(i + 1   ) * GRID_COLS + (j - 1  )] - old_val;
            }

            ap_fixed<8, 1> D;
            if(i == GRID_ROWS - 1                          ) {
                D = 0;
            }
            else {
                D = imgvf[(i + 1   ) * GRID_COLS + (j      )] - old_val;
            }

            ap_fixed<8, 1> DR;
            if(i == GRID_ROWS - 1  ||  j == GRID_COLS - 1  ) {
                DR = 0;
            }
            else {
                DR = imgvf[(i + 1   ) * GRID_COLS + (j + 1  )] - old_val;
            }

            ap_fixed<8, 1> a = heaviside(UL);
            ap_fixed<8, 1> b = heaviside(U);
            ap_fixed<8, 1> c = heaviside(UR);
            ap_fixed<8, 1> d = heaviside(L);
            ap_fixed<8, 1> e = heaviside(R);
            ap_fixed<8, 1> f = heaviside(DL);
            ap_fixed<8, 1> g = heaviside(D);
            ap_fixed<8, 1> h = heaviside(DR);



            ap_fixed<8, 1> tmp = MU_O_LAMBDA;

            ap_fixed<8, 1> vHe = old_val + tmp * (a * UL + b * U + c * UR + d * L + e * R + f * DL + g * D + h * DR);

            ap_fixed<8, 1> vI = I[i * GRID_COLS + j];

            ap_fixed<8, 1> tmp2 = ONE_O_LAMBDA;

            ap_fixed<8, 1> new_val = vHe - (tmp2 * vI * (vHe - vI));
            
            result[i * GRID_COLS + j] = new_val;
    
            
            ap_fixed<8, 1> absolute = new_val - old_val;
            if(absolute < 0) {
                absolute = -absolute;
            }

            total_diff += absolute;
        }
    }

    //return 1;

    return (total_diff / (ap_fixed<8, 1>)(GRID_ROWS * GRID_COLS)).to_float(); // makes it crash
}

void workload(ap_fixed<8, 1> result[GRID_ROWS * GRID_COLS], ap_fixed<8, 1> imgvf[GRID_ROWS * GRID_COLS], ap_fixed<8, 1> I[GRID_ROWS * GRID_COLS])
{
    #pragma HLS INTERFACE m_axi port=result offset=slave bundle=result1
    #pragma HLS INTERFACE m_axi port=imgvf offset=slave bundle=imgvf1
    #pragma HLS INTERFACE m_axi port=I offset=slave bundle=I1
    
    #pragma HLS INTERFACE s_axilite port=result bundle=control
    #pragma HLS INTERFACE s_axilite port=imgvf bundle=control
    #pragma HLS INTERFACE s_axilite port=I bundle=control
    
    #pragma HLS INTERFACE s_axilite port=return bundle=control


    int i;
    float diff = 1.0;

    for (i = 0; i < ITERATION / 2; i++) {
        diff = lc_mgvf(result, imgvf, I);
        diff = lc_mgvf(imgvf, result, I);
    }
    return;

}
}
