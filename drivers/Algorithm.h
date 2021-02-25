/** \file Algorithm.h
 *
 *  Contains the declarations of functions used in Algorithm.c file
 */

#ifndef Algorithm_h
#define Algorithm_h

#include <stdio.h>

extern "C" {

    /**
 * Executes 1st step of the Sensor Fusion Algorithm.
 * Produces a 1D array which is the Support Degree Matrix when given
 * a 1D array of sensor readings at one time stamp and the number
 * of sensors being considered.
 */
double* sdm_calculator(double[],int);

/**
 * Executes a part of 2nd step of the Sensor Fusion Algorithm.
 * Produces a 1D array consisting of EigenValues for the given
 * support degree matrix in the form of a 1D array and size of that array
 */
double* eigen_value_calculation(double*,int);

/**
 * Executes a part of 2nd step of the Sensor Fusion Algorithm.
 * Produces a 1D array consisting of EigenVectors of one
 * EigenValue when support degree matrix in 1D array form,
 * the array size and the ordinal position of the EigenValue
 * when arranged in descending order is given to the function.
 */
double* eigen_vector_calculation(double*,int,int);

/**
 * Computes contribution rate of each sensor at a specific
 * time stamp when supplied with EigenValues in descending order
 * and the number of sensors in consideration.
 */
double* compute_alpha(double[],int);

/**
 * Computes accumulated contribution rates of sensors at a specific
 * time stamp when supplied with contribution rates and the number
 * of sensors in consideration.
 */
double* compute_phi(double[], int);

/**
 * Calculates the integrated support degree score of each sensor
 */
double* compute_integrated_support_degree_score(double[], double[], double[], double[],double, int);

/**
 * Calculates the fused value after identifying and removing
 * faulty sensor values.
 */
double faulty_sensor_and_sensor_fusion(double[],double[],double, int);

}


#endif /* SensorFusionAlgorithm_h */
