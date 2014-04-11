/*
 *  a_model.h
 *  diffuser_z
 *
 *  Created by Arkadiy Simonov on 8/12/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

extern double two_times(double);
extern void multDoubleArray(double *x,int size);
extern void calculate_scattering(double *scattering,double* inp_params);


#ifdef __cplusplus
} /* end extern "C" */
#endif


