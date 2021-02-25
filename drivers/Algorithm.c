/** \file SensorFusionAlgorithm.c
 *  Contains functions performing each step of the algorithm.
 *  The calculation of Eigen Values, Eigen Vectors and integrated support
 *  degree score has been implemented using an external library called gsl.
 */

#include "Algorithm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

/** \brief Calculate Support Degree Matrix.
 *
 *  Creates the support degree matrix as a 1D array instead of a 2D array
 *  to make it easier for the calculation of Eigen Values and Eigen Vectors
 *  using the gsl library.
 *
 *  @param[in] sensorinputs Readings of all sensors for a specific timestamp.
 *  @param[in] size The number of sensors being considered.
 *
 *  \return A 1D array which is the Support Degree Matrix.
 */
double* sdm_calculator(double sensorinputs[], int size){

    int i,j,k=0;
    //printf("size %d\n", size);
    double *dmatrix = (double *) malloc(sizeof(double)*(size*size));

    for(i=0; i<size;i++){

        for(j=0;j<size;j++){
            double temp = sensorinputs[i]-sensorinputs[j];

            if(temp<0){
                temp = -temp;
            }

            temp = exp(-(temp));
            dmatrix[k]= temp;
            //printf("%lf\t",dmatrix[k]);
            k++;

            if(k%size==0){
            //printf("\n");
            }
        }
    }
    return dmatrix;
}

/** \brief Calculates EigenValues for a given Support Degree Matrix.
 *
 *  Creates an 1D array of EigenValues for a given support degree matrix
 *  using functions from the gsl external library.
 *
 *  @param[in] dmatrix The Support Degree matrix produced by sdm_calculator
 *   function.
 *  @param[in] size The number of sensors being considered which is same as
 *   size of dmatrix.
 *
 *  \return A 1D array of EigenValues in descending order for the given dmatrix.
 */
double* eigen_value_calculation(double* dmatrix, int size){
    int i;
    gsl_matrix_view m = gsl_matrix_view_array (dmatrix, size, size);
    gsl_vector *eval = gsl_vector_alloc (size);
    gsl_matrix *evec = gsl_matrix_alloc (size, size);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (size);
    gsl_eigen_symmv (&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
    double *eval_i = (double *) malloc(sizeof(double)*(size));
    for(i = 0; i < size; i++){
        eval_i[i] = gsl_vector_get (eval, i);
        //printf("eigenvalue = %g\n", eval_i[i]);
    }
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    return eval_i;
}

/** \brief Produces EigenVectors for a specific EigenValue of
 *   Support Degree Matrix.
 *
 *  Creates an 1D array of EigenVectors for a given support degree matrix
 *  using functions from the gsl external library. The actual EigenValue is
 *  not supplied but the position of EigenValue in a descending sorted
 *  array of EigenValues as returned by eigen_value_calculation method is
 *  given.
 *
 *  @param[in] dmatrix The Support Degree matrix produced by sdm_calculator
 *   function.
 *  @param[in] size The number of sensors being considered which is same as
 *   size of dmatrix.
 *  @param[in] column The column number of the EigenVector matrix produced
 *   by the gsl library. This is same as the position of EigenValue when
 *   the EigenValues are sorted in descending order.
 *
 *  \return A 1D array of EigenVectors for the given EigenValue.
 */
double* eigen_vector_calculation(double* dmatrix, int size, int column){
    int i;
    gsl_matrix_view m = gsl_matrix_view_array (dmatrix, size, size);
    gsl_vector *eval = gsl_vector_alloc (size);
    gsl_matrix *evec = gsl_matrix_alloc (size, size);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (size);
    gsl_eigen_symmv (&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
    double *evec_i =(double *) malloc(sizeof(double)*(size));
    for(i=0;i<size;i++){
        evec_i[i]= gsl_matrix_get(evec, i, column);
        //printf("Eigen Vector: %g\n", evec_i[i]);
    }
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    return evec_i;
}

/** \brief Calculates the contribution rates i.e. alphas of every sensor
 *  for a specific time stamp.
 *
 * Produces an array of contribution rates a.k.a alphas where the contribution
 * rate of nth sensor is the nth element in that array.
 *
 * @param[in] eval_i EigenValues in descending order as produced by
 *  eigen_value_calculation function.
 * @param[in] size Number of EigenValues which is same as the
 *   number of sensors and number of contribution rate values.
 *
 * \return Contribution rates of all sensors for a specific time stamp.
 */
double* compute_alpha(double eval_i[],int size){
    double *list_of_alphas = (double *) malloc(sizeof(double)*(size));
    double sum_of_evals = 0.0;

    //Summation of EigenValues
    for(int i=0; i<size; i++){
        sum_of_evals += eval_i[i];
    }
    //Constructing an array of alphas a.k.a contribution rates
    for(int i=0;i <size;i++){
        list_of_alphas[i] = eval_i[i]/sum_of_evals;
        //printf("list of alpha = %lf\n\n",list_of_alphas[i]*100);
    }
    return list_of_alphas;
}

/** \brief Calculates the accumulated contribution rates i.e. phis for a specific
 *   time stamp.
 *
 * Produces an array of accumulated contribution rates where the accumulated
 * contribution rate of first n sensors is the nth element in that array.
 *
 * @param[in] list_of_alphas Contribution rates i.e., alphas as produced by
 *  compute_alpha function.
 * @param[in] size Number of contribution rates which is the same as the
 *   number of sensors.
 *
 * \return Array of accumulated contribution rates of n sensors
 *   for specific time stamp.
 */
double* compute_phi(double list_of_alphas[], int size){
    double *list_of_phi = (double *) malloc(sizeof(double)*(size));

    //The Contribution rate and the Accumulated Contribution rate of first
    //sensor is the same
    list_of_phi[0] = list_of_alphas[0];
    //printf("list of phi = %lf\n",list_of_phi[0]*100);

    for(int i=1; i<size;i++){
        list_of_phi[i] = list_of_phi[i-1] + list_of_alphas[i];
        //printf("list of phi = %lf\n",list_of_phi[i]*100);
    }
    return list_of_phi;
}

/** \brief Calculates the integrated support degree score of all sensors at a
 * specific time stamp.
 *
 * Produces an array of integrated support degree scores where the score
 * nth sensor is the nth element in that array.
 *
 * @param[in] sensorinputs Readings of all sensors at a specific timestamp.
 * @param[in] list_of_alphas Contribution rates a.k.a alphas where the
 *  contribution rate of nth sensor is the nth element in that array.
 * @param[in] list_of_phi Accumulated contribution rates a.k.a phis where
 *  the accumulated contribution rate of m sensors is the mth element in array.
 * @param[in] dmatrix Support Degree Matrix of the sensor readings at a
 *  specific time stamp.
 * @param[in] criterion The minimum value of accumulated contribution rate
 *  required to determine the number of principal components to be taken
 *  into consideration. The expected value is not in percentage but between
 *  0 and 1.
 * @param[in] size The number of sensors.
 *
 * \return An array consisting of the integrated support degree score
 *   of each sensor at a specific timestamp.
 */
double* compute_integrated_support_degree_score(double sensorinputs[],
            double list_of_alphas[], double list_of_phi[], double dmatrix[],
            double criterion, int size){

    double calculation,temp=0,temp1=0,temp2=0;
    gsl_matrix *y_temp =  gsl_matrix_alloc(size, size);
    gsl_matrix *y =  gsl_matrix_alloc(size, size);
    double *list_of_m_phi = (double *) malloc(sizeof(double)*(size));
    double *evec_i = (double *) malloc(sizeof(double)*(size));
    int i,j,k,l=0,rows,col,o,p,q,rows1=0,col1=0,rows2,col2;
    for(i=0;i<size;i++){
        list_of_m_phi[i]=list_of_phi[i];
        if(list_of_phi[i]>criterion){
            break;
        }
    }
    gsl_matrix *dmatrix2d =  gsl_matrix_alloc(size, size);

    for (j = 0; j < size; j++){
        for (k = 0; k < size; k++){
            gsl_matrix_set(dmatrix2d, j, k, dmatrix[l++]);
            //printf ("m(%d,%d) = %g\n", j, k, gsl_matrix_get (dmatrix2d, j, k));
        }
    }

    for(o=0;o<size;o++){
        evec_i = eigen_vector_calculation(sdm_calculator(sensorinputs,size),size,o);
        for(col=0;col<size;col++){
            for(rows=0;rows<size;rows++){
                calculation = evec_i[col]* gsl_matrix_get(dmatrix2d,rows,col);
                gsl_matrix_set(y_temp, rows, col, calculation);
                //printf ("m(%d,%d) = %g\n", rows, col, gsl_matrix_get (y_temp, rows, col));
            }
        }
        for(rows=0;rows<size;rows++){
            temp = 0;
            for(col=0;col<size;col++){
                temp += gsl_matrix_get (y_temp, rows, col);
                //printf("TEMP:%lf \n",temp);
            }

        if(rows1==size){
            col1++;
            rows1=0;
        }

        gsl_matrix_set(y, rows1++, col1, temp);
        //printf("TEMP cols:%lf \n",temp);
        }
    }

    for(col=0;col<size;col++){
        for(rows=0;rows<size;rows++){
            //printf ("m(%d,%d) = %g\n", rows, col, gsl_matrix_get (y, rows, col));
        }
    }
    gsl_matrix *Z_temp =  gsl_matrix_alloc(size, size);
    double *Z = (double *) malloc(sizeof(double)*(size));
    for(p=0;p<size;p++){
        for(q=0;q<size;q++){
            temp1 =list_of_alphas[p]* gsl_matrix_get(y,q,p);
            gsl_matrix_set(Z_temp, q, p, temp1);
            //printf ("m(%d,%d) = %g\n", p, q, gsl_matrix_get (Z_temp, q, p));
            //printf("alpha :%lf \n",list_of_alphas[q]);
        }
    }
    for(rows2=0;rows2<size;rows2++){
        temp2 = 0;
        for(col2=0;col2<=i;col2++){
        temp2 += gsl_matrix_get (Z_temp, rows2, col2);
        //printf("TEMP:%lf \n",temp2);
        }

        if(rows2==size){
            col2++;
            rows2=0;
        }
        Z[rows2] = temp2;
    }
 //   for(rows2=0;rows2<size;rows2++){
  //      printf ("m(%d) = %lf\n", rows2, Z[rows2] );
  //  }
  gsl_matrix_free(y_temp);
  gsl_matrix_free(y);
  gsl_matrix_free(dmatrix2d);
  gsl_matrix_free(Z_temp);
  free(list_of_m_phi);
  free(evec_i);
  return Z;
}

/** \brief Determines a fused reading by eliminating erroneous readings
 *
 * Outputs a list of faulty sensors and a final fused reading value
 * for a specific time stamp on to a text file.
 *
 * @param[in] Z Array containing integrated support degree score of each
 *  sensor at a specific timestamp.
 * @param[in] inputsensors Readings of all sensors for a specific timestamp.
 * @param[in] criterion from the user multiplying to the average.
 * @param[in] size The number of sensors being considered.
 *
 * \return The fused reading value after eliminating faulty sensor readings.
 */
double faulty_sensor_and_sensor_fusion(double Z[], double inputsensors[], double criterion,
		int size){

    int i, tempfault=0,j=0;
    double *weight = (double *) malloc(sizeof(double)*(size));
    int *fault = (int *) malloc(sizeof(int)*(size));
    double average, sum=0,calculation=0,fusion_value=0;

    //Summation of integrated support degree score of all sensors
    for(i=0;i<size;i++){
        sum += Z[i];
    }

    average = fabs((sum/size))*criterion;
   // printf("Criterion: %lf",criterion);
    //printf("Average : %lf", average);
    for(i=0;i<size;i++){

    	//Identifying a faulty sensor and storing it's index in an int array called fault
        if(fabs(Z[i])<average){
            tempfault = i;
          //  printf("Fault Detected! The sensor number %d is a faulty sensor!\n",tempfault+1);
            fault[j]=tempfault;
            j++;
        }
	}

	//Making the reading and score of all faulty sensors to zero
    for(i=0;i<j;i++){
        Z[fault[i]]=0;
        inputsensors[fault[i]]=0;
    }
   // for(i=0;i<size;i++){
   //     printf("Z = %lf\n",Z[i]);
   // }

    //Assigning the sum of integrated support degree score of all
    //NON-FAULTY sensors to 'calculation' variable
    for(i=0;i<size;i++){
        calculation += Z[i];
    }

    //Creating an array of weight coefficients of each sensor by dividing a
    //sensor's integrated support degree score by sum of all scores
    for(i=0;i<size;i++){
        weight[i] = Z[i]/calculation;
    }
  //  for(i=0;i<size;i++){
  //      printf("Weight coefficient : %lf\n",weight[i]);
  //  }

    //Calculating the fused value as a summation of product of
    //weight coefficient and sensor reading
    for(i=0;i<size;i++){
        fusion_value += weight[i] * inputsensors[i];
    }
    free(weight);
    free(fault);
    return fusion_value;
}
