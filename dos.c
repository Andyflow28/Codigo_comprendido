#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Qromb.h"


float A[100], B[100];
double DOS,f,g, normalizacion;
double En = -0.4;
double Hop_t = -0.2;
/*
double En= 0.4;
double Hop_t=-0.2;
*/
int N=100;
int n;
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif
double Func1(double x);
double Func2(double x);
double dos_f(double x);



/***********************************************************************************************************************/

int main(void)  {

    double LowerLimitDx = M_PI- acos(En/(4*Hop_t));
    double Delta0 = 24.0*sqrt(2.0);

    FILE *fp;

    char* file_names[] = {"c0g00/c0g001.dat", "c0g00/c0g005.dat", "c0g00/c0g010.dat", "c0g00/c0g015.dat", "c0g00/c0g020.dat"};
    char* DOS_names[] = {"c0g_output/c0g001_output.dat", "c0g_output/c0g005_output.dat", "c0g_output/c0g010_output.dat", "c0g_output/c0g015_output.dat", "c0g_output/c0g020_output.dat"};

    for(int i=0; i<5; i++){
        fp = fopen(file_names[i], "r");
        for(int j=0 ; j<N ; j++){
            fscanf(fp, "%f", &A[j]);            
            fscanf(fp, "%f", &B[j]);
        }
        fclose(fp);
        
        FILE *dos;
        printf("- Se ha creado con exito el archivo: %s\n",DOS_names[i]);
        dos=fopen(DOS_names[i],"w");
        for(n=0; n<N; n++){
            normalizacion = qromb(dos_f, LowerLimitDx, M_PI);					
            f= qromb(Func1,LowerLimitDx, M_PI);
            g= qromb(Func2,LowerLimitDx, M_PI);
            DOS = (f + g)/(normalizacion);
            fprintf(dos,"%f\t %f\n", A[n], DOS);
        }
        fclose(dos);
        if (dos == NULL) {
            printf("--- Error en %s\n", DOS_names[i]);
            return 1;
        }
        printf("+ Se han guardados los datos correctamente en: %s\n\n", DOS_names[i]);
    }
}
/********************************************************************************************************************/
double Func1(double x)
{
    double  result, aux3, aux4, rho_k, a_k, b, positive_k, Delta0 = 24.0*sqrt(2.0);

    float aux1, z, aux2, vx, vz, cos_x, sin_x, jacobian;
    float vel, cos_z, sin_z, der_zx, factor;

    cos_x = cos(x);
    sin_x = sqrt(1-cos_x*cos_x);
    aux1 = En + 2*Hop_t*cos_x;

    aux2 = Hop_t;
    z = M_PI-acos(0.5*aux1/aux2);
    cos_z = cos(z);
    sin_z = sqrt(1-cos_z*cos_z);
    vx = 2*Hop_t*sin_x;
    vz = 2*Hop_t*sin_z;
    vel = sqrt(vx*vx + vz*vz);
    der_zx = pow(vx/vz,2);
    factor = sqrt(1 + der_zx);
    jacobian = factor/vel;      

/*
    a_k = pow(A[n],2.0) - pow(B[n],2.0) - pow(Delta0*(cos_x-cos_z),2.0);
*/
    a_k = pow(A[n],2.0) - pow(B[n],2.0) - pow(cos_x-cos_z,2.0);
    b = 2*A[n]*B[n];
    rho_k = sqrt( pow(a_k,2.0) + pow(b,2.0) );
    positive_k = 1 + a_k/rho_k;
    aux3 = sqrt(2.0*rho_k);
    aux4 = A[n]/aux3;	
    result = jacobian*aux4*sqrt(positive_k);
    return result;
}

/***********************************************************************************************************************/
double Func2(double x)
{	
    double  result, aux3, aux4, rho_k, a_k, b, negative_k, Delta0 = 24.0*sqrt(2.0);

    float aux1, z, aux2, vx, vz, cos_x, sin_x, jacobian;
    float vel, cos_z, sin_z, der_zx, factor;

    cos_x = cos(x);
    sin_x = sqrt(1-cos_x*cos_x);
    aux1 = En + 2*Hop_t*cos_x;

    aux2 = Hop_t;
    z = M_PI-acos(0.5*aux1/aux2);
    cos_z = cos(z);
    sin_z = sqrt(1-cos_z*cos_z);
    vx = 2*Hop_t*sin_x;
    vz = 2*Hop_t*sin_z;
    vel = sqrt(vx*vx + vz*vz);
    der_zx = pow(vx/vz,2);
    factor = sqrt(1 + der_zx);
    jacobian = factor/vel;      


/*
    a_k = pow(A[n],2.0) - pow(B[n],2.0) - pow(Delta0*(cos_x-cos_z),2.0);
*/
    a_k = pow(A[n],2.0) - pow(B[n],2.0) - pow(cos_x-cos_z,2.0);
    b = 2*A[n]*B[n];
    rho_k = sqrt( pow(a_k,2.0) + pow(b,2.0) );
    negative_k = 1 - a_k/rho_k;
    aux3 = sqrt(2.0*rho_k);
    aux4 = B[n]/aux3;	
    result = jacobian*aux4*sqrt(negative_k);
    return result;
}
/***********************************************************************************************************************/
double dos_f(double x)  
{
    float aux1, z, aux2, vx, vz, cos_x, sin_x, jacobian;
    float vel, cos_z, sin_z, der_zx, factor;

    cos_x = cos(x);
    sin_x = sqrt(1-cos_x*cos_x);
    aux1 = En + 2*Hop_t*cos_x;

    aux2 = Hop_t;	
    z = M_PI-acos(0.5*aux1/aux2);
    cos_z = cos(z);
    sin_z = sqrt(1-cos_z*cos_z);
    vx = 2*Hop_t*sin_x;
    vz = 2*Hop_t*sin_z;
    vel = sqrt(vx*vx + vz*vz);
    der_zx = pow(vx/vz,2);
    factor = sqrt(1 + der_zx);
    jacobian = factor/vel;      
    return jacobian;
}

