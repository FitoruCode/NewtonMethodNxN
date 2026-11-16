#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "nrutil.h"

typedef struct newton_D
{
	int n;      // dimension
	int Nmax;   // max iterations
	int k;      // iteration number
	double eps; // precision stop criterion
	double delta; // numerical derivative aproximation step'
	double gamma; // step scaler TODO - alg for it
	double yNorm; // previous step norm
	double currentNorm; // current step norm
	
	double *x0; // R^n x0 vector
	double *y0; // R^n f(x0) vector
	double *(*f)(double*, double*, int); // R^n->R^n function
	double **D; // Derivative matrix
	double **Dinv; // Derivative matrix inverse
	double *dk; // step vectors
	double *hv; // h search vector

	FILE *fp; // output file
}NEWTON_D;

int main(void);
double *f(double *x, double *y, int n);
double *set_x0(double *x, int n);
void set_NEWTON_D(NEWTON_D *s);
void free_memory(NEWTON_D *s);
void print_dvector(double *x, int n);
void print_dmatrix(double **m, int n);
double **get_derivative_matrix_value(double *(*f)(double *x, double *y, int n), double *x, double *y, double **m, double delta, int n);
double **inv(double **m, double **inv, int n);
double *multiply_matrix_v(double **m, double *x, double *y, int n);
double *scale_vector(double *x, double gamma, int n);
double vector_norm(double *x, int n);
double *subtract_vector(double *x, double *y, int n);
bool stop_criterion(double normPrevious, double normCurrent, double eps, int Nmax, int k);
void print_step(double *x, double yNorm, int k, int n);
void log_step(double *x, double yNorm, int k, int n, FILE* fp);
int search_j(double *(*f)(double *x, double *y, int n), double *x, double *y, double *dk, double *hv, double C, int n);
double matrix_norm(double **x, int n);
int classical_newton(void);
double get_h_constant(double **D, double *d, int n);


int main()
{
	classical_newton();
	return 0;
}

int classical_newton(void)
{
	NEWTON_D s;
	set_NEWTON_D(&s);

	while(true)
	{
		print_step(s.x0, s.yNorm, s.k, s.n);
		log_step(s.x0, s.yNorm, s.k, s.n, s.fp);
		
		s.D = get_derivative_matrix_value(s.f, s.x0, s.y0, s.D, s.delta, s.n);
		s.Dinv = inv(s.D, s.Dinv, s.n);
		s.dk = multiply_matrix_v(s.Dinv, s.y0, s.dk, s.n);
		s.dk = scale_vector(s.dk, s.gamma, s.n);
		s.x0 = subtract_vector(s.x0, s.dk, s.n);
		s.y0 = f(s.x0, s.y0, s.n);
		s.currentNorm = vector_norm(s.y0, s.n);
		
		if (stop_criterion(s.yNorm, s.currentNorm, s.eps, s.Nmax, s.k)) {break;}
		s.yNorm = s.currentNorm;
		s.k += 1;
	}

	free_memory(&s);
	return 0;
}

// Example function from end of chapter 12 in lab pdf script
double *f(double *x, double *y, int n)
{
	y[0] =  pow((pow(x[0], 2.0)+pow(x[1], 2.0)), 2.0) - 2*((pow(x[0], 2.0)-pow(x[1], 2.0)));
	y[1] =  pow((pow(x[0], 3.0)+pow(x[1], 3.0)), 2.0) + pow(x[1], 3.0) - pow(x[0]-0.5, 2.0); 
	return y;
}

int search_j(double *(*f)(double *x, double *y, int n), double *x, double *y, double *dk, double *hv, double C, int n)
{
	double tau;
	double h0=0.0;
	double htau;
	int j=1;

	for(int8_t i=0; i<n; i++)
	{
		h0 += y[i]*y[i];
	}

	while(true)
	{
		tau = pow(2, -j);
		for(int8_t i=0; i<n; i++)
		{
			hv[i] = x[i]-tau*dk[i];
		}
		y=f(hv,y,n);

		htau = 0.0;
		for(int8_t i=0; i<n; i++)
		{
			htau += y[i]*y[i];
		}
		if(htau <= h0-pow(2, -j)*C)
		{
			return j;
		}
		j=j=+1;
	}
	return -1;
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
	s->gamma=0.1;
	s->x0=dvector(0, s->n-1);
	s->y0=dvector(0, s->n-1);
	s->D=dmatrix(0, s->n-1, 0, s->n-1);
	s->Dinv=dmatrix(0, s->n-1, 0, s->n-1);
	s->dk=dvector(0, s->n-1);
	s->hv=dvector(0, s->n-1);
	
	s->Nmax=100;
	s->eps=1e-6;
	s->delta=1e-5;

	s->x0=set_x0(s->x0, s->n);
	s->y0=s->f(s->x0, s->y0,s->n);
	s->yNorm=vector_norm(s->y0, s->n);
	s->currentNorm=s->yNorm+99;

	s->fp=fopen("./output/out.txt", "w");
}

void free_memory(NEWTON_D *s)
{
	free_dvector(s->x0, 0, s->n-1);
	free_dvector(s->y0, 0, s->n-1);
	free_dvector(s->dk, 0, s->n-1);
	free_dvector(s->hv, 0, s->n-1);
	free_dmatrix(s->D, 0, s->n-1, 0, s->n-1);
	free_dmatrix(s->Dinv, 0, s->n-1, 0, s->n-1);
}

bool stop_criterion(double normPrevious, double normCurrent, double eps, int Nmax, int k)
{
	if (k>=Nmax) {return true;}
	if (fabs(normPrevious-normCurrent)<eps) {return true;}
	return false;
}

double **get_derivative_matrix_value(double *(*f)(double *x, double *y, int n), double *x, double *y, double **m, double delta, int n)
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
		x[i] += delta;
		y = f(x,y,n);
		for(int8_t j=0; j<n; j++)
		{
			m[j][i] += y[j];
			m[j][i] /= delta;
		}
		x[i] -= delta;
	}
	return m;
}

double *subtract_vector(double *x, double *y, int n)
{
	for(int8_t i=0; i<n; i++)
	{
		x[i] -= y[i];
	}
	return x;
}

double vector_norm(double *x, int n)
{
	double norm = 0.0;
	for(int8_t i=0; i<n; i++)
	{
		norm += x[i]*x[i];
	}
	return sqrt(norm);
}

double matrix_norm(double **x, int n)
{
	double norm = 0.0;
	for(int8_t i=0; i<n; i++)
	{
		for(int8_t j=0; j<n; j++)
		{
			norm += x[i][j]*x[i][j];
		}
	}
	return sqrt(norm);
}

double get_h_constant(double **D, double *d, int n)
{
	double result=matrix_norm(D, n)*vector_norm(d, n)*0.25;
	return result;
}

double *scale_vector(double *x, double gamma, int n)
{
	for(int8_t i=0; i<n; i++)
	{
		x[i] = x[i]*gamma;
	}
	return x;
}

void log_step(double *x, double yNorm, int k, int n, FILE* fp)
{
	if (fp == NULL)
    {
        return;
    }
	fprintf(fp, "Step %d: x={", k);
	for(int8_t i=0; i<n-1; i++)
	{
		fprintf(fp, "%f, ", x[i]);
	}
	fprintf(fp, "%f}", x[n-1]);
	fprintf(fp, ", ||y||=%f", yNorm);
	fprintf(fp, "\n");
}

void print_step(double *x, double yNorm, int k, int n)
{
	printf("Step %d: x={", k);
	for(int8_t i=0; i<n-1; i++)
	{
		printf("%f, ", x[i]);
	}
	printf("%f}", x[n-1]);
	printf(", ||y||=%f", yNorm);
	printf("\n");
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

double *multiply_matrix_v(double **m, double *x, double *y, int n)
{
	double v;
	for(int8_t i=0; i<n; i++)
	{
		v = 0.0;
		for(int8_t j=0; j<n; j++)
		{
			v = v+m[i][j]*x[j];
		}
		y[i] = v;
	}
	return y;
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