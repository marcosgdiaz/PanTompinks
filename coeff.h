#ifndef COEFF_H
#define COEFF_H

#include <arm_math.h>
/* Generated with MatlabSOS2CMSIS
Phil Birch, University of Sussex, 2017*/

#define NUM_SECTIONS_BANDPASS 7
#define NUM_SECTIONS_DERIVATIVE 1
#define NUM_SECTIONS_MOVING 15

#endif
/*Example usage:
#include "coeff"
float32_t pState[NUM_SECTIONS*4]={0};
arm_biquad_casd_df1_inst_f32 S;

In the your main function init the filter
arm_biquad_cascade_df1_init_f32(&S,NUM_SECTIONS,pCoeffs,pState);
To run the filter:
arm_biquad_cascade_df1_f32(&S,pSrc,pDest,BUFFER_SIZE);
See CMSIS doc for varible descriptions
*/
