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
	double **D; // Derivative matrix
	double **Dinv; // Derivative matrix inverse

}NEWTON_D;

int main(void);
double *f(double *x, double *y, int n);
double *set_x0(double *x, int n);
void set_NEWTON_D(NEWTON_D *s);
void print_dvector(double *x, int n);
void print_dmatrix(double **m, int n);
double **get_derivative_matrix_value(double *(*f)(double *x, double *y, int n), double *x, double *y, double **m, double theta, int n);
double **inv(double **m, double **inv, int n);

int main(void)
{
	NEWTON_D s;
	set_NEWTON_D(&s);

	print_dvector(s.x0, s.n);
	print_dvector(s.y0, s.n);

	get_derivative_matrix_value(s.f, s.x0, s.y0, s.D, s.delta, s.n);
	inv(s.D, s.Dinv, s.n);

	print_dmatrix(s.D, s.n);
	print_dmatrix(s.Dinv, s.n);

	return 0;
}

// Example function from end of chapter 12 in lab pdf script
double *f(double *x, double *y, int n)
{
	y[0] =  pow((pow(x[0], 2.0)+pow(x[1], 2.0)), 2.0) - 2*((pow(x[0], 2.0)-pow(x[1], 2.0)));
	y[1] =  pow((pow(x[0], 3.0)+pow(x[1], 3.0)), 2.0) + pow(x[1], 3.0) - pow(x[0]-0.5, 2.0); 
	return y;
}

// Example x0
double *set_x0(double *x, int n)
{
	x[0] = 1;
	x[1] = 1;
	return x;
}

void set_NEWTON_D(NEWTON_D *s)
{
	s->k=0;
	s->f=f;
	s->n=2;
	s->x0=dvector(0, s->n);
	s->y0=dvector(0, s->n);
	s->D=dmatrix(0, s->n, 0, s->n);
	s->Dinv=dmatrix(0, s->n, 0, s->n);
	
	s->Nmax=10;
	s->eps=1e-6;
	s->delta=1e-5;

	s->x0=set_x0(s->x0, s->n);
	s->y0=s->f(s->x0, s->y0,s->n);
}

double **get_derivative_matrix_value(double *(*f)(double *x, double *y, int n), double *x, double *y, double **m, double theta, int n)
{
	for(int8_t i=0; i<n; i++)
	{
		for(int8_t j=0; j<n; j++)
		{
			m[i][j] = -y[i];
		}
	}

	for(int8_t i=0; i<n; i++)
	{
		x[i] += theta;
		y = f(x,y,n);
		for(int8_t j=0; j<n; j++)
		{
			m[j][i] += y[j];
			m[j][i] /= theta;
		}
		x[i] -= theta;
	}
	return m;
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

void print_dmatrix(double **m, int n)
{
	for(int8_t i=0; i<n; i++)
	{
		print_dvector(m[i], n);
	}
}

double **inv(double **m, double **inv, int n)
{
    if (n <= 0 || m == NULL || inv == NULL)
        return NULL;
    double **a = malloc(n * sizeof(double *));

    if (a == NULL)
        return NULL;

    for (int i = 0; i < n; i++)
    {
        a[i] = malloc(n * sizeof(double));

        if (a[i] == NULL)
        {
            for (int k = 0; k < i; k++)
                free(a[k]);
            free(a);
            return NULL;
        }
        for (int j = 0; j < n; j++)
            a[i][j] = m[i][j];
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = (i == j) ? 1.0 : 0.0;

    for (int i = 0; i < n; i++)
    {
        int pivot = i;

        for (int j = i + 1; j < n; j++)
            if (fabs(a[j][i]) > fabs(a[pivot][i]))
                pivot = j;

        if (fabs(a[pivot][i]) < 1e-12)
        {
            for (int k = 0; k < n; k++)
                free(a[k]);

            free(a);
            return NULL;
        }

        for (int j = 0; j < n; j++)
        {
            double temp = a[i][j];
            a[i][j] = a[pivot][j];
            a[pivot][j] = temp;

            temp = inv[i][j];
            inv[i][j] = inv[pivot][j];
            inv[pivot][j] = temp;
        }

        double div = a[i][i];
        for (int j = 0; j < n; j++)
        {
            a[i][j] /= div;
            inv[i][j] /= div;
        }

        for (int j = 0; j < n; j++)
        {
            if (j == i)
                continue;

            double factor = a[j][i];
            for (int k = 0; k < n; k++)
            {
                a[j][k] -= factor * a[i][k];
                inv[j][k] -= factor * inv[i][k];
            }
        }
    }

    for (int i = 0; i < n; i++)
        free(a[i]);
    free(a);
    return inv;
}