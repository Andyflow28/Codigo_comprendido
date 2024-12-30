// Qromb.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Qromb.h"
#include "Utils.h"


double qromb(double (*func)(double x), double a, double b);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double trapzd(double (*func)(double x), double a, double b, int n);

double qromb1(double (*func)(double x), double a, double b);
void polint1(double xa[], double ya[], int n, double x, double *y, double *dy);
double trapzd1(double (*func)(double x), double a, double b, int n);

double qromb2(double (*func)(double x), double a, double b);
void polint2(double xa[], double ya[], int n, double x, double *y, double *dy);
double trapzd2(double (*func)(double x), double a, double b, int n);

double qromb3(double (*func)(double x), double a, double b);
void polint3(double xa[], double ya[], int n, double x, double *y, double *dy);
double trapzd3(double (*func)(double x), double a, double b, int n);

double qromb4(double (*func)(double x), double a, double b);
void polint4(double xa[], double ya[], int n, double x, double *y, double *dy);
double trapzd4(double (*func)(double x), double a, double b, int n);

double qromb5(double (*func)(double x), double a, double b);
void polint5(double xa[], double ya[], int n, double x, double *y, double *dy);
double trapzd5(double (*func)(double x), double a, double b, int n);

static double xsav;
static double (*nrfunc)(double,double);


//************************************************************************************

#define EPS 1.0E-7
#define JMAX 50
#define JMAXP (JMAX+1)
#define K 5



double qromb(double (*func)(double x), double a, double b)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) {
			// Invoke the next line to print the number of iterations.
			// printf("interations = %d\n",j); 
					 return ss;}
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}

double qromb1(double (*func)(double x), double a, double b)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd1(func,a,b,j);
		if (j >= K) {
			polint1(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) {
			// Invoke the next line to print the number of iterations.
			// printf("interations = %d\n",j); 
					 return ss;}
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb1");
	return 0.0;
}

double qromb2(double (*func)(double x), double a, double b)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd2(func,a,b,j);
		if (j >= K) {
			polint2(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) {
			// Invoke the next line to print the number of iterations.
			// printf("interations = %d\n",j); 
					 return ss;}
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb2");
	return 0.0;
}
double qromb3(double (*func)(double x), double a, double b)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd3(func,a,b,j);
		if (j >= K) {
			polint3(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) {
			// Invoke the next line to print the number of iterations.
			// printf("interations = %d\n",j); 
					 return ss;}
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb3");
	return 0.0;
}
double qromb4(double (*func)(double x), double a, double b)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd4(func,a,b,j);
		if (j >= K) {
			polint4(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) {
			// Invoke the next line to print the number of iterations.
			// printf("interations = %d\n",j); 
					 return ss;}
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb4");
	return 0.0;
}
double qromb5(double (*func)(double x), double a, double b)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd5(func,a,b,j);
		if (j >= K) {
			polint5(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) {
			// Invoke the next line to print the number of iterations.
			// printf("interations = %d\n",j); 
					 return ss;}
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb5");
	return 0.0;
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K


//************************************************************************************


void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}

void polint1(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint1");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}

void polint2(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint2");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}

void polint3(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint3");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}

void polint4(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint4");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}

void polint5(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint5");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}


//************************************************************************************


#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double x), double a, double b, int n)
{
        double x,tnm,sum,del;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}

double trapzd1(double (*func)(double x), double a, double b, int n)
{
        double x,tnm,sum,del;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}

double trapzd2(double (*func)(double x), double a, double b, int n)
{
        double x,tnm,sum,del;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}

double trapzd3(double (*func)(double x), double a, double b, int n)
{
        double x,tnm,sum,del;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}

double trapzd4(double (*func)(double x), double a, double b, int n)
{
        double x,tnm,sum,del;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}

double trapzd5(double (*func)(double x), double a, double b, int n)
{
        double x,tnm,sum,del;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}

#undef FUNC
