#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Qromb.h"


float A[100], B[100];
double DOS, plus, minus, term, normalizacion;
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
double Termino_positivo(double x);
double Termino_negativo(double x);
double calcular_jacobiano(double x);
double dos_f(double x);
double Calculo_de_los_g0();
double Calculo_de_los_g0_sin_desorden();
double Calculo_dat();
double Calculo_data_sin_desorden();

/***********************************************************************************************************************/

int main(void)  {
    //!Calculo_de_los_g0();
    Calculo_dat();
    //!Calculo_data_sin_desorden();
    //!Calculo_de_los_g0_sin_desorden();
}
/********************************************************************************************************************///? Calculo con desorden
double Termino_positivo(double x)
{
    double result, rho_k, a_k, b_k, positive_k, Delta0 = 24.0*sqrt(2.0);
    float cos_x, z, cos_z, jacobian;

    // Calcular el jacobiano 
    jacobian = calcular_jacobiano(x);

    cos_x = cos(x);
    z = M_PI - acos(0.5 * (En + 2 * Hop_t * cos_x) / Hop_t);
    cos_z = cos(z);

    // Precalcular los términos a_k y b_k 
    a_k = (A[n] * A[n]) - (B[n] * B[n]) - (cos_x - cos_z) * (cos_x - cos_z);
    b_k = 2 * A[n] * B[n];
    rho_k = sqrt(a_k * a_k + b_k * b_k);
    positive_k = 1 + a_k / rho_k;

    result = jacobian * (A[n] / sqrt(2.0 * rho_k)) * sqrt(positive_k);
    return result;
}

/************************************************************* ********************///! Calculo del Jacobiano
double calcular_jacobiano(double x)
{
    float aux1, z, vx, vz, cos_x, sin_x;
    float vel, cos_z, sin_z, der_zx, factor;

    cos_x = cos(x);
    sin_x = sqrt(1 - cos_x * cos_x);
    aux1 = En + 2 * Hop_t * cos_x;

    z = M_PI - acos(0.5 * aux1 / Hop_t);
    cos_z = cos(z);
    sin_z = sqrt(1 - cos_z * cos_z);
    
    vx = 2 * Hop_t * sin_x;
    vz = 2 * Hop_t * sin_z;
    vel = sqrt(vx * vx + vz * vz);
    der_zx = (vx / vz) * (vx / vz);  
    factor = sqrt(1 + der_zx);

    return factor / vel;
}

/*************************************************************** *///! 

double Termino_negativo(double x)
{	
    double result, rho_k, a_k, b_k, negative_k, Delta0 = 24.0*sqrt(2.0);
    float cos_x, z, cos_z, jacobian;

    jacobian = calcular_jacobiano(x);

    cos_x = cos(x);
    z = M_PI - acos(0.5 * (En + 2 * Hop_t * cos_x) / Hop_t);
    cos_z = cos(z);

    // Precalcular los términos a_k y b_k de la forma más eficiente
    a_k = (A[n] * A[n]) - (B[n] * B[n]) - (cos_x - cos_z) * (cos_x - cos_z);
    b_k = 2 * A[n] * B[n];
    rho_k = sqrt(a_k * a_k + b_k * b_k);
    negative_k = 1 - a_k / rho_k;

    result = jacobian * (B[n] / sqrt(2.0 * rho_k)) * sqrt(negative_k);
    return result;
}

/***********************************************************************************************************************/
double dos_f(double x)  
{
    double jacobian;

    // Calcular el jacobiano 
    jacobian = calcular_jacobiano(x);

    return jacobian;
}

//? //////////////////////////////////////////////////////////////////////////////////////////////////////

double Calculo_de_los_g0() {
    double LowerLimitDx = M_PI- acos(En/(4*Hop_t));
    double Delta0 = 24.0*sqrt(2.0);

    FILE *fp;

    char* file_names[] = {"c0g00/c0g001.dat", "c0g00/c0g005.dat", "c0g00/c0g010.dat", "c0g00/c0g015.dat", "c0g00/c0g020.dat"};
    char* DOS_names[] = {"c0g_output/c0g001_output.dat", "c0g_output/c0g005_output.dat", "c0g_output/c0g010_output.dat", "c0g_output/c0g015_output.dat", "c0g_output/c0g020_output.dat"};

    for(int i=0; i<5; i++){
        fp = fopen(file_names[i], "r");

        if (fp == NULL) {
        printf("Error al abrir el archivo: %s\n", file_names[i]);
        return 1;
        }

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
            plus= qromb(Termino_positivo,LowerLimitDx, M_PI);
            minus= qromb(Termino_negativo,LowerLimitDx, M_PI);
            DOS = (plus + minus)/(normalizacion);
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

//todo //////////////////////////////////////////////////////////////////////////////////////////// Seccion sin desorden

double Termino_Sin_desorden(double x) {
    double result, rho_k, a_k, Delta0 = 24.0*sqrt(2.0);
    float omega_N, k_xz, term_sin_desorden, lower_k, cos_x, z, cos_z, factor, jacobian;

    // Calcular el jacobiano usando la nueva función
    jacobian = calcular_jacobiano(x);

    cos_x = cos(x);
    z = M_PI - acos(0.5 * (En + 2 * Hop_t * cos_x) / Hop_t);
    cos_z = cos(z);

    k_xz = cos_x - cos_z;
    omega_N = A[n];
    lower_k = omega_N * omega_N - k_xz * k_xz;
    term_sin_desorden = omega_N / sqrt(lower_k);

    result = jacobian * term_sin_desorden;
    return result;
}

double Calculo_de_los_g0_sin_desorden() {
    double LowerLimitDx = M_PI- acos(En/(4*Hop_t));
    double Delta0 = 24.0*sqrt(2.0);

    FILE *fp;

    char* file_names[] = {"c0g00/c0g001.dat", "c0g00/c0g005.dat", "c0g00/c0g010.dat", "c0g00/c0g015.dat", "c0g00/c0g020.dat"};
    char* DOS_names[] = {"data_sin_desorde/c0g001_output_sin_desorde.dat", "data_sin_desorde/c0g005_output_sin_desorde.dat", "data_sin_desorde/c0g010_output_sin_desorde.dat", "data_sin_desorde/c0g015_output_sin_desorde.dat", "data_sin_desorde/c0g020_output_sin_desorde.dat"};

    for(int i=0; i<5; i++){
        fp = fopen(file_names[i], "r");

        if (fp == NULL) {
        printf("Error al abrir el archivo: %s\n", file_names[i]);
        return 1;
        }

        for(int j=0 ; j<N ; j++){
            fscanf(fp, "%f", &A[j]);            
            fscanf(fp, "%*f");
        }
        fclose(fp);
        
        FILE *dos;
        printf("- Se ha creado con exito el archivo: %s\n",DOS_names[i]);
        dos=fopen(DOS_names[i],"w");
        for(n=0; n<N; n++){
            normalizacion = qromb(dos_f, LowerLimitDx, M_PI);					
            term= qromb(Termino_Sin_desorden,LowerLimitDx, M_PI);
            DOS = term/(normalizacion);
            // for (int i = 0; i < n; i++) {
            //     printf("%f\n", A[i]);
            // }
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

//! ///////////////////// Sistema de calculo con los datos no G ????????????????????????///////////////

double Calculo_dat() {
    double LowerLimitDx = M_PI- acos(En/(4*Hop_t));
    double Delta0 = 24.0*sqrt(2.0);

    FILE *fp;

    char* file_names[] = {"data/c00g000E04t02.dat", "data/c00g001E04t02.dat", "data/c00g005E04t02.dat", "data/c00g010E04t02.dat", "data/c00g015E04t02.dat"};
    char* DOS_names[] = {"output/output_000E04t02.dat", "output/output_001E04t02.dat", "output/output_005E04t02.dat", "output/output_010E04t02.dat", "output/output_015E04t02.dat"};

    for(int i=0; i<5; i++){
        fp = fopen(file_names[i], "r");

        if (fp == NULL) {
        printf("Error al abrir el archivo: %s\n", file_names[i]);
        return 1;
        }

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
            plus= qromb(Termino_positivo,LowerLimitDx, M_PI);
            minus= qromb(Termino_negativo,LowerLimitDx, M_PI);
            DOS = (plus + minus)/(normalizacion);
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

//todo /////////////////////////////////// el mismo calculo de arriba pero sin el desorden //////////


double Calculo_data_sin_desorden() {
    double LowerLimitDx = M_PI- acos(En/(4*Hop_t));
    double Delta0 = 24.0*sqrt(2.0);

    FILE *fp;

    char* file_names[] = {"data/c00g000E04t02.dat", "data/c00g001E04t02.dat", "data/c00g005E04t02.dat", "data/c00g010E04t02.dat", "data/c00g015E04t02.dat"};
    char* DOS_names[] = {"output_sin_desorden/output_sin_desorden_000E04t02.dat", "output_sin_desorden/output_sin_desorden_001E04t02.dat", "output_sin_desorden/output_sin_desorden_005E04t02.dat", "output_sin_desorden/output_sin_desorden_010E04t02.dat", "output_sin_desorden/output_sin_desorden_015E04t02.dat"};

    for(int i=0; i<5; i++){
        fp = fopen(file_names[i], "r");

        if (fp == NULL) {
        printf("Error al abrir el archivo: %s\n", file_names[i]);
        return 1;
        }

        for(int j=0 ; j<N ; j++){
            fscanf(fp, "%f", &A[j]);            
            fscanf(fp, "%*f");
        }
        fclose(fp);
        
        FILE *dos;
        printf("- Se ha creado con exito el archivo: %s\n",DOS_names[i]);
        dos=fopen(DOS_names[i],"w");
        for(n=0; n<N; n++){
            normalizacion = qromb(dos_f, LowerLimitDx, M_PI);					
            term= qromb(Termino_Sin_desorden,LowerLimitDx, M_PI);
            DOS = term/(normalizacion);
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