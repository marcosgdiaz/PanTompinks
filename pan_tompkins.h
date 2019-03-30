/*
 * thresholding.h
 *
 *  Created on: Feb 22, 2019
 *      Author: marcosgonzalezdiaz
 */

#ifndef SRC_PAN_TOMPKINS_H_
#define SRC_PAN_TOMPKINS_H_

#define SAMPLE_LENGTH 2000 // Testing sample length
#define BLOCK_SIZE 100 //Processing splitted in blocks
#define F_SAMPLING 200 //Frequency sampling
#define MAX_PEAKS 25 //Maximum number of peaks during a block size
#define DELAY 30 //0.15*F_SAMPLING Twice the delay between filtered and final signal approximately
#define SLOPE 15 //0.075*F_SAMPLING
#define QRS 8 //maximum number of qrs stored
#define MISSED 40 //0.2*F_SAMPLING from Matlab code

#include <arm_math.h> //define ARM_MATH_M4 and add corresponding libraries and paths
#include <stdio.h>
#include <stdlib.h>
#include <nrf_delay.h>
#include <math.h>
#include <SEGGER_RTT.h>
#include "coeff.h" //where coefficients of the filters are located

void FindPeaks(float32_t* HR_final, uint32_t minimumDistance, uint32_t* number_peaks, uint32_t index);
void SignalProcessing(float32_t* HR, uint32_t i);
void Filters_init(void);
void Threshold_init(void);
void Thresholding(uint32_t index);
void LocatePeak(uint32_t loc, uint32_t index, bool possible);
void Diff(uint32_t* vector, uint32_t size, float32_t* out);
void Append(uint32_t* locs, float32_t* peaks, uint32_t loc, float32_t peak, uint32_t* n_pks);
float32_t Abs(float32_t num);
void SaveBuffer(void);
void RateCalculator(void);

#endif /* SRC_PAN_TOMPKINS_H_ */
