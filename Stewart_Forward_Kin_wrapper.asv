
/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#define SIMPLIFIED_RTWTYPES_COMPATIBILITY
#include "rtwtypes.h"
#undef SIMPLIFIED_RTWTYPES_COMPATIBILITY
#endif



/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
# ifndef MATLAB_MEX_FILE
# include <Arduino.h>
// change
int i,j,k,l,b;
float p=0,q=0,sumf=0,pin=3.141592,tolf=1e-7,tola=1e-7,maxiter=100;
int num;
double a[6] = {0,0,400,0,0,0};
double topoffset[3] = {0,0,62};
double baseoffset[3] = {0,0,57};
double func[6],rotationmatrix[3][3],A[6][6],B[6],platform[6],L[6][6],U[6][6],Y[6],deltaa[6],al[6];
float array[3],array1[3],array2[3],array3[3],bas[3],xB[3];



# endif
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 6
#define u_1_width 4
#define y_width 6

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output function
 *
 */
void Stewart_Forward_Kin_Outputs_wrapper(const real_T *actual_lengths,
			const real_T *platform_parameters,
			real_T *actual_pos,
			const real_T *xD)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
# ifndef MATLAB_MEX_FILE
// change

while (num<maxiter) {
        
    for (i=0;i<6;i++){
        al[i]=actual_lengths[i]+389.9;
    }

    float sumf=0,p=0,q=0;
    
    double gamma=a[3]*pin/180;
    double beta=a[4]*pin/180;
    double alpha=a[5]*pin/180;
    
    rotationmatrix[0][0] = cos(gamma)*cos(beta);
    rotationmatrix[0][1] = cos(gamma)*sin(beta)*sin(alpha)-sin(gamma)*cos(alpha);
    rotationmatrix[0][2] = cos(gamma)*sin(beta)*cos(alpha)+sin(gamma)*sin(alpha);
    rotationmatrix[1][0] = sin(gamma)*cos(beta);
    rotationmatrix[1][1] = sin(gamma)*sin(beta)*sin(alpha)+cos(gamma)*cos(alpha);
    rotationmatrix[1][2] = sin(gamma)*sin(beta)*cos(alpha)-cos(gamma)*sin(alpha);
    rotationmatrix[2][0] = -sin(beta);
    rotationmatrix[2][1] = cos(beta)*sin(alpha);
    rotationmatrix[2][2] = cos(beta)*cos(alpha);
     
     
    array2[0] = rotationmatrix[0][0]*topoffset[0]+rotationmatrix[0][1]*topoffset[1]+rotationmatrix[0][2]*topoffset[2];
    array2[1] = rotationmatrix[1][0]*topoffset[0]+rotationmatrix[1][1]*topoffset[1]+rotationmatrix[1][2]*topoffset[2];
    array2[2] = rotationmatrix[2][0]*topoffset[0]+rotationmatrix[2][1]*topoffset[1]+rotationmatrix[2][2]*topoffset[2];

    array3[0]=array2[0]+baseoffset[0];
    array3[1]=array2[1]+baseoffset[1];
    array3[2]=array2[2]+baseoffset[2];
    
    // adjust x, y, z
    double x=a[0];
    double y=a[1];
    double z=a[2];
    x=x-array3[0];
    y=y-array3[1];
    z=z-array3[2];

    for (i=0;i<6;i++){
        b=i%2;
        if (b==0) {
            p=(60*(i+1))-(platform_parameters[1]/2);
            q=(60*(i+1))-(platform_parameters[2]/2);
        }else{
            p=platform_parameters[1]+p;
            q=platform_parameters[2]+q;
        }
        
        
        bas[0]=platform_parameters[3]*cos((q/180.0)*pin);
        bas[1]=platform_parameters[3]*sin((q/180.0)*pin);
        bas[2]=0;
        
        xB[0]=x-bas[0];
        xB[1]=y-bas[1];
        xB[2]=z-bas[2];

        
        array[0]=platform_parameters[0]*cos((p/180.0)*pin);
        array[1]=platform_parameters[0]*sin((p/180.0)*pin);
        array[2]=0;
        
        array1[0] = rotationmatrix[0][0]*array[0]+rotationmatrix[0][1]*array[1]+rotationmatrix[0][2]*array[2];
        array1[1] = rotationmatrix[1][0]*array[0]+rotationmatrix[1][1]*array[1]+rotationmatrix[1][2]*array[2];
        array1[2] = rotationmatrix[2][0]*array[0]+rotationmatrix[2][1]*array[1]+rotationmatrix[2][2]*array[2];

        // compute fi(a)
        func[i] = pow((xB[0]+array1[0]),2) + pow((xB[1]+array1[1]),2) + pow((xB[2]+array1[2]),2) - pow(al[i],2);
        sumf = sumf + fabs(func[i]); // add to sum of fi(a)
        B[i] = -func[i];
        
        // compute A partial derivatives
        A[i][0]=2*(xB[0]+array1[0]);
        A[i][1]=2*(xB[1]+array1[1]);
        A[i][2]=2*(xB[2]+array1[2]);
        A[i][3]=2*(-xB[0]*array1[1]+xB[1]*array1[0]);
        A[i][4]=2*(array1[2]*(-xB[0]*cos(gamma)+xB[1]*sin(gamma))-xB[2]*(array[0]*cos(beta)+array[1]*sin(beta)*sin(alpha)));
        A[i][5]=2*array[1]*(xB[0]*rotationmatrix[0][2]+xB[1]*rotationmatrix[1][2]+xB[2]*rotationmatrix[2][2]);
    }
    
    if (sumf < tolf){
        // printf("sumf");
        break;
    }

    // LU decomp step, solve for delta a (j): A*deltaa=B
    // double* deltaa = LUD(A,B);

    for (l = 0; l < 6; l++) {
        // Upper Triangular
        for (k = l; k < 6; k++) {
            // Summation of L(i, j) * U(j, k)
            double sum = 0.0;
            for (j = 0; j < l; j++)
                sum += (L[l][j] * U[j][k]);
            // Evaluate U(i, k)
            U[l][k] = A[l][k] - sum;
        }
        // Lower Triangular
        for (k = l; k < 6; k++) {
            if (l == k)
                L[l][l] = 1; // Diagonal as 1
            else {
                // Summation of L(k, j) * U(j, i)
                double sum = 0.0;
                for (j = 0; j < l; j++)
                    sum += (L[k][j] * U[j][l]);
                // Evaluate L(k, i)
                L[k][l] = (A[k][l] - sum) / U[l][l];
            }
        }
    }

    for(l=0; l<6; l++)
    {
        Y[l]=B[l];
        for(j=0; j<l; j++)
        {
            Y[l]-=L[l][j]*Y[j];
        }
    }

    for(l=5; l>=0; l--)
    {
        deltaa[l]= Y[l];
        for(j=l+1; j<6; j++)
        {
            deltaa[l]-=U[l][j]*deltaa[j];
        }
        deltaa[l]/=U[l][l];
    }
    // deltaa[0]=0.01;
    // deltaa[1]=0.01;
    // deltaa[2]=0.01;
    // deltaa[3]=0.01;
    // deltaa[4]=0.01;
    // deltaa[5]=0.01;

    deltaa[3]=deltaa[3]*180/pin;
    deltaa[4]=deltaa[4]*180/pin;
    deltaa[5]=deltaa[5]*180/pin; // convert back to deg
    double sumd=0;
    for (i=0;i<6;i++){
        sumd+=fabs(deltaa[i]);
    }

    if (sumd < tola){
        // printf("sumd");
        break;
    }
    
    
    
    num+=1;
    for (i=0;i<6;i++){
        a[i]=a[i]+deltaa[i];
    }
}
    
    float temp=a[3];
    a[3] = a[5];
    a[5]=temp;
    for (i=0;i<6;i++){
        actual_pos[i]=a[i];
    }



    return actual_pos;

#endif
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}

/*
 * Updates function
 *
 */
void Stewart_Forward_Kin_Update_wrapper(const real_T *actual_lengths,
			const real_T *platform_parameters,
			real_T *actual_pos,
			real_T *xD)
{
/* %%%-SFUNWIZ_wrapper_Update_Changes_BEGIN --- EDIT HERE TO _END */
 
/* %%%-SFUNWIZ_wrapper_Update_Changes_END --- EDIT HERE TO _BEGIN */
}

