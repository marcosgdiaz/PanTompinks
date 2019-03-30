/*
 * thresholding.c
 *
 *  Created on: Feb 22, 2019
 *      Author: marcosgonzalezdiaz
 */

#include "pan_tompkins.h"

//Filtering variables
float32_t HR_filtered[BLOCK_SIZE];
float32_t HR_derivative[BLOCK_SIZE]; //It can be removed for saving memory
float32_t HR_squared[BLOCK_SIZE];  //It can be removed for saving memory
float32_t HR_final[BLOCK_SIZE];
float32_t HR_buffer[DELAY]; //Saving the last samples of the last processed block
static float32_t pState_bandpass[NUM_SECTIONS_BANDPASS*4]={0};  //NUM_SECTIONS is defined in coeff.h
static float32_t pState_derivative[NUM_SECTIONS_DERIVATIVE*4]={0};  //NUM_SECTIONS is defined in coeff.h
static float32_t pState_moving[NUM_SECTIONS_MOVING*4]={0};  //NUM_SECTIONS is defined in coeff.h
static arm_biquad_casd_df1_inst_f32 Bandpass, Derivative, Moving;   //structure that defines the filter
extern float32_t pCoeffs_bandpass[],pCoeffs_derivative[],pCoeffs_moving[]; //coefficients of the filters
//*************

//Thresholding variables
static bool thr_fetched = false;
static bool thr_fetched_1 = false;
static float32_t thr_sig, sig_level, thr_sig_1, sig_level_1, thr_noise, noise_level, thr_noise_1, noise_level_1;
static float32_t m_selected_RR = 0;
uint32_t number_peaks = 0;
static uint32_t number_qrs = 0;
uint32_t number_qrs_1 = 0;
static uint32_t qrs_loc[QRS]={0};
uint32_t qrs_loc_1[QRS]={0};
static float32_t qrs_peak[QRS]={0};
float32_t qrs_peak_1[QRS]={0};
uint8_t HeartRate = 0;
//*************

//Variables related to peaks
static float32_t memory = 0;
static bool possible = false;
uint32_t locs[MAX_PEAKS] = {0};
float32_t peaks[MAX_PEAKS] = {0};
uint32_t locs_1[MAX_PEAKS] = {0};
float32_t peaks_1[MAX_PEAKS] = {0};
//*******

/**
  @brief Both dynamic thresholds are initialized
  @note It can be improved since the BLOCK_SIZE may be small for the initialization
  */

void Threshold_init(void){
	uint32_t aux;

	if(thr_fetched==false){
		arm_max_f32(HR_final,BLOCK_SIZE,&thr_sig,&aux);
		thr_sig *= 0.3333;
		sig_level = thr_sig;
		arm_mean_f32(HR_final,BLOCK_SIZE,&thr_noise);
		thr_noise *= 0.5;
		noise_level = thr_noise;
		thr_fetched = true;
	}

	if(thr_fetched_1==false){
		arm_max_f32(HR_filtered,BLOCK_SIZE,&thr_sig_1,&aux);
		thr_sig_1 *= 0.3333;
		sig_level_1 = thr_sig_1;
		arm_mean_f32(HR_filtered,BLOCK_SIZE,&thr_noise_1);
		thr_noise_1 *= 0.5;
		noise_level_1 = thr_noise_1;
		thr_fetched_1 = true;
	}
}

/**
  @brief It initializes all the filters
  @note There are three filters: Bandpass filter, Derivative filter and
  Moving-window integration filter
  */

void Filters_init(void){
	arm_biquad_cascade_df1_init_f32(&Bandpass,NUM_SECTIONS_BANDPASS,pCoeffs_bandpass,pState_bandpass);
	arm_biquad_cascade_df1_init_f32(&Derivative,NUM_SECTIONS_DERIVATIVE,pCoeffs_derivative,pState_derivative);
	arm_biquad_cascade_df1_init_f32(&Moving,NUM_SECTIONS_MOVING,pCoeffs_moving,pState_moving);
}
//*******

/**
  @brief It processes the blocks in order to calculate the Heart Rate
  @param  HR    the array that contains the raw waveform. It must be length BLOCK_SIZE
  @param  i     i-th block that is being processed
  @note you can find what each step is doing in the function definition
  */

void SignalProcessing(float32_t* HR, uint32_t i){
	//Bandpass Filter
	arm_biquad_cascade_df1_f32(&Bandpass, HR, HR_filtered,BLOCK_SIZE);
	arm_scale_f32(HR_filtered, 0.0000597, HR_filtered, BLOCK_SIZE); //Amplitude gain of the filter
	//Derivative Filter
	arm_biquad_cascade_df1_f32(&Derivative, HR_filtered, HR_derivative, BLOCK_SIZE);
	//Square
	arm_mult_f32(HR_derivative, HR_derivative, HR_squared, BLOCK_SIZE);
	//Moving Window
	arm_biquad_cascade_df1_f32(&Moving, HR_squared, HR_final, BLOCK_SIZE);
	arm_scale_f32(HR_final, 0.033, HR_final, BLOCK_SIZE); //Amplitude gain of the filter
	if(i==1){
		Threshold_init();
	}
	if(thr_fetched && thr_fetched_1){
		//Seeking all possible QRS complex
		FindPeaks(HR_final, MISSED, &number_peaks,i);
		//Determining whether the peaks are QRS complex or not
		Thresholding(i);
		//number_peaks = 0; Discomment when the visualitation ends
	}
	//Saving the last 30 samples of the filtered signal
	SaveBuffer();
	//Calculating the Heart Rate
	RateCalculator();
}

/**
  @brief It finds peaks and select the maximum one if they are within a minimumDistance
  @param  HR_final	the array that contains the signal after all filters
  @param  minimumDistance	the minimum distance that two peaks can be found
  @param number_peaks	the address od the variable that contains the number of peaks that have been found
  @param  index		i-th block that is being processed
  @note number_peaks must be reset to zero after the thresholding
  */

void FindPeaks(float32_t* HR_final, uint32_t minimumDistance, uint32_t* number_peaks, uint32_t index){
	for(uint32_t i = 0; i < BLOCK_SIZE ; i++){

		//The first sample of the block may be a peak
		if(possible){
			if(memory > HR_final[i]){
				locs[*(number_peaks)] = BLOCK_SIZE*index + i;//A -1 is missing on purpose
				peaks[*(number_peaks)] = HR_final[i];
				LocatePeak(locs[*(number_peaks)], *(number_peaks), possible);
				*(number_peaks)+=1;
			}
			possible = false;
		}

		//Determine whether the sample is a peak or not
		if(memory < HR_final[i] && HR_final[i+1] < HR_final[i] && index + i > 0 && i != BLOCK_SIZE - 1){
			//If the distance between peaks is less that the minimum distance, the maximum peak will be selected
			if(*(number_peaks) > 0 && BLOCK_SIZE*index + i - locs[*(number_peaks)-1] < minimumDistance &&
					HR_final[i] > peaks[*(number_peaks)-1]){
				locs[*(number_peaks)-1] = BLOCK_SIZE*index + i;
				peaks[*(number_peaks)-1] = HR_final[i];
				LocatePeak(locs[*(number_peaks)-1], *(number_peaks)-1,possible);
			}
			//Otherwise the peak is added directly
			else if((*(number_peaks) == 0 || BLOCK_SIZE*index + i - locs[*(number_peaks)-1] > minimumDistance)
					&& *(number_peaks) < MAX_PEAKS){
				locs[*(number_peaks)] = BLOCK_SIZE*index + i;
				peaks[*(number_peaks)] = HR_final[i];
				LocatePeak(locs[*(number_peaks)], *(number_peaks), possible);
				*(number_peaks)+=1;
			}
		}
		//Saving the last sample for using it in the next block
		else if(i == BLOCK_SIZE - 1 && HR_final[i] > memory){
			possible = true;
		}
		memory = HR_final[i];
	}
}

/**
  @brief Determine whether the peak corresponds to a QRS complex or not and update the thresholds
  @param index	i-th block that is being processed
  @note the search back is not completely implemented and a more simple solution has been done
  @note T-wave detection has been removed due to its low importance (from my point of view)
  */

void Thresholding(uint32_t index){
	for(uint32_t i=0; i < number_peaks; i++){

		//Calculate the mean of the RR intervals
		if(number_qrs == QRS){
			float32_t diffRR[QRS-1];
			Diff(qrs_loc, QRS, diffRR);
			arm_mean_f32(diffRR,QRS-1,&m_selected_RR);
			float32_t interval = qrs_loc[QRS] - qrs_loc[QRS-1];
			if(interval <= 0.92*m_selected_RR || interval >= 1.16*m_selected_RR){
				thr_sig *= 0.5;
				thr_sig_1 *= 0.5;
			}
		}

		//Determine whether a QRS complex is missing or not
		if(m_selected_RR!=0){
			if(locs[i] - qrs_loc[number_qrs-1] >= 1.66*m_selected_RR){
				//Intead of seaching back a peak (complex task),
				//a peak is assumed to be at m_selected_RR distance from the last QRS complex
				Append(qrs_loc, qrs_peak, qrs_loc[number_qrs-1] + m_selected_RR, 0, &number_qrs);
				Append(qrs_loc_1, qrs_peak_1, qrs_loc_1[number_qrs-1] + m_selected_RR, 0, &number_qrs_1);
			}
		}

		//Determine whether the peak is QRS complex or not
		if(peaks[i] >= thr_sig){
			Append(qrs_loc, qrs_peak, locs[i], peaks[i], &number_qrs);
			if(peaks_1[i] >= thr_sig_1){
				Append(qrs_loc_1, qrs_peak_1,locs_1[i], peaks_1[i], &number_qrs_1);
				sig_level_1 = 0.125 * peaks_1[i] + 0.875 * sig_level_1;
			}
			sig_level = 0.125 * peaks[i] + 0.875 * sig_level;
		}
		else{
			noise_level = 0.125 * peaks[i] + 0.875 * noise_level;
			noise_level_1 = 0.125 * peaks_1[i] + 0.875 * noise_level_1;
		}

		//Uploading dynamic thresholds
		if(noise_level != 0 || sig_level != 0){
			thr_sig = noise_level + 0.25*(Abs(sig_level - noise_level));
			thr_noise = 0.5 * thr_sig;
		}
		if(noise_level_1 != 0 || sig_level_1 != 0){
			thr_sig_1 = noise_level_1 + 0.25*(Abs(sig_level_1 - noise_level_1));
			thr_noise_1 = 0.5 * thr_sig_1;
		}
	}
}

/**
  @brief Locate certain peak of the final signal in the filtered signal
  @param  loc	location of the peak of the final signal that has to be located in the filtered signal
  @param  index	number of peaks found at that moment
  @param  possible whether the index location of the peak is zero (true) or not (false)
  @note Depending of the variable loc, it might access to the buffer
  */

void LocatePeak(uint32_t loc, uint32_t index, bool possible){
	float32_t aux,aux_1;
	uint32_t l1, l2;
	int32_t comp = loc%BLOCK_SIZE - DELAY;
	//Locate in the current block
	if(comp > 0){
		arm_max_f32(&HR_filtered[loc%BLOCK_SIZE - DELAY], DELAY, &aux,(locs_1+index));
		*(peaks_1+index) = aux;
		*(locs_1+index) += loc - DELAY;//adding the initial delay
	}
	//Locate in both current block and buffered block
	else{
		arm_max_f32(HR_filtered, loc%BLOCK_SIZE+1, &aux, &l1);
		arm_max_f32(&HR_buffer[loc%BLOCK_SIZE], DELAY - loc%BLOCK_SIZE, &aux_1, &l2);
		if(aux > aux_1){
			*(peaks_1+index) = aux;
			*(locs_1+index) = loc - loc%BLOCK_SIZE + l1;
		}
		else if(possible){
			*(peaks_1+index) = aux_1;
			*(locs_1+index) = loc - DELAY + l2 - 1;
		}
		else{
			*(peaks_1+index) = aux_1;
			*(locs_1+index) = loc - DELAY + l2;
		}
	}
}

/**
  @brief calculate the difference vector defined as out: out[i] = in[i+1] - in[i]
  @param  vector	input vector
  @param  size		size of the input vector
  @param  out	output vector
  @note it is used to calculate the RR intervals
  */

void Diff(uint32_t* vector, uint32_t size, float32_t* out){
	for(uint32_t i=0; i < size - 1; i++){
		*(out+i) = *(vector+i+1) - *(vector+i);
	}
}

/**
  @brief First Input First Output in order to save the last 8 QRS complex
  @param  locs	qrs peak location vector
  @param  peaks	qrs peak vector
  @param  loc qrs peak location that will be added in the vector
  @param  peak qrs peak that will be added in the vector
  @param  n_pks number of qrs complex at that moment (maximum 8)
  @note If the distance between two consecutive qrs complex is less than 40 samples
   then the vector may be updated depending on the peak value
  */

void Append(uint32_t* locs, float32_t* peaks, uint32_t loc, float32_t peak, uint32_t* n_pks){
	//Case the distance between two qrs is less than MISSED. Uploading with the maximum peak
	if(loc - *(locs + *n_pks - 1) < MISSED && *n_pks > 0 && peak > *(peaks + *n_pks - 1)){
		*(locs + *n_pks - 1) = loc;
		*(peaks + *n_pks - 1) = peak;
	}

	else if(*n_pks == QRS && loc - *(locs + *n_pks - 1) > MISSED){
		for(uint32_t i=0; i < QRS - 1; i++){
			*(locs+i) = *(locs + i + 1);
			*(peaks+i) = *(peaks + i + 1);
		}
		*(peaks + QRS - 1) = peak;
		*(locs + QRS - 1) = loc;
	}

	else if(*n_pks==0 || loc - *(locs + *n_pks - 1) > MISSED){
		*(peaks + *n_pks) = peak;
		*(locs + *n_pks) = loc;
		*(n_pks)+=1;
	}
}

/**
  @brief Saving last 30 samples in order to be able to locate peaks in the filtered signal
  @note Only the filtered signal is buffered
  */

void SaveBuffer(void){
	for(uint32_t i = 0; i<DELAY; i++){
		HR_buffer[i] = HR_filtered[BLOCK_SIZE-1-DELAY+i];
	}
}

/**
  @brief Absolute value of a float32_t number
  */

float32_t Abs(float32_t num){
	if(num < 0) return num*-1;
	else return num;
}

/**
  @brief Calculates the heart rate knowing the last 8 QRS complex
  @note It is not robust to non-wanted spikes derived from the sensor, so it can be improved
  */

void RateCalculator(void){
	if(number_qrs_1 == QRS){
		float32_t diffRR[QRS-1];
		float32_t aux;
		Diff(qrs_loc_1, QRS, diffRR);
		arm_mean_f32(diffRR,QRS-1,&aux);
		HeartRate = (uint8_t) (F_SAMPLING*60/aux);
	}
}
