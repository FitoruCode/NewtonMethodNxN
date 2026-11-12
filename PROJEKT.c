#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

typedef struct newton_D
{
	int n;      // dimension
	int Nmax;   // max iterations
	int k;      // iteration number
	double eps; // precision stop criterion
	double delta; // numerical derivative aproximation step
	
	double *x0; // R^n x0 vector
	double *y0; // R^n f(x0) vector
	double *(*f)(double*, double*, int); // R^n->R^n function

}NEWTON_D;

int main(void);
double *f(double *x, double *y, int n);
double *set_x0(double *x, int n);
void set_NEWTON_D(NEWTON_D *s);
void print_dvector(double *x, int n);

int main(void)
{
	NEWTON_D s;
	set_NEWTON_D(&s);

	print_dvector(s.x0, s.n);
	print_dvector(s.y0, s.n);

	return 0;
}

// Example function from end of chapter 12 in lab pdf script
double *f(double *x, double *y, int n)
{
	y[0] =  pow((pow(x[0], 2.0)+pow(x[1], 2.0)), 2.0) - 2*((pow(x[0], 2.0)-pow(x[1], 2.0)));
	y[1] =  pow((pow(x[0], 3.0)+pow(x[1], 3.0)), 3.0) + pow(x[1], 3.0) - pow(x[0]-0.5, 2.0); 
	return y;
}

// Example x0
double *set_x0(double *x, int n)
{
	x[0] = 1;
	x[1] = 2;
	return x;
}

void set_NEWTON_D(NEWTON_D *s)
{
	s->k=0;
	s->f=f;
	s->x0=dvector(0, s->n);
	s->y0=dvector(0, s->n);
	
	s->n=2;
	s->Nmax=10;
	s->eps=1e-6;
	s->delta=1e-5;

	s->x0=set_x0(s->x0, s->n);
	s->y0=s->f(s->x0, s->y0,s->n);
}

void print_dvector(double *x, int n)
{
	for (int8_t i=0; i<n; i++)
	{
		printf("%f", x[i]);
		printf(" ");
	}
	printf("\n");
}