#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define pi 3.14

#define _Mi    40      //[g/mol]
#define _Z     1  
#define _Ei    50      //[eV]
#define _dEi   0*_Ei //[eV]
#define _alpha 0*pi  //[rad]

#define _h    0.010  //[m]  
#define _dh   0.001  //[m]
#define _l    0.000 //[m]

#define _H    0.014      //[m]
#define _L    _l + 0.030 //[m]
#define _dL   0.001      //[m]

#define _mi (_Mi*0.001/6e23) //[kg]
#define _qi (_Z*1.6e-19)     //[Kl]

#define _dt   1e-8 //[s]
#define _N    1e5  //[ions]
#define _Umax 100  //[V]

#define _fileName "50eV_slit00mm.txt"

double zero_one_dist();
double uniform(double start, double stop);
double gauss(double E, double Var);

void slit_start(double *t, double *rV);
void slit_f    (double *t, double *rV, double *F);
void slit_step (double *t, double *rV);
int  slit_check(double *rV);
int  slit_track(double *t, double *rV);

void track_f    (double U, double *t, double *rV, double *F);
void track_step (double U, double *t, double *rV);
int  track_check(double *rV);
int  track_track(double U, double *t, double *rV);

int main()
{
    srand(clock());

    FILE *file;
    file = fopen(_fileName, "w");

    for (double U = 0; U < _Umax; U += 1)
    {
        printf("U = %.1f\n", U);
        
        double t, rV[4];
        int count = 0;
        for (int i = 0; i < _N; i++)
        {
            if (slit_track(&t, rV) == 1)
            {
                if (track_track(U, &t, rV) == 1)
                    count++;
            }
        }
        fprintf(file, "%.1f\t%d\t%f\n", U, count, float(count/_N));
    }

    fclose(file);

    printf("\nHello world\n");
    return 0;
}

//===========================
//====    Runge-Kutta    ====
//===========================

//=============================
//====    Entrance slit    ====
//=============================
void slit_start(double *t, double *rV)
{
    *t  = 0;
    *(rV + 0) = 0;
    *(rV + 1) = uniform(_h - _dh, _h + _dh);
    double a = uniform(0.5*(pi - _alpha), 0.5*(pi + _alpha));
    double Vi = sqrt(2*gauss(_Ei, _dEi)*1.6e-19/_mi);
    *(rV + 2) = Vi*sin(a);
    *(rV + 3) = Vi*cos(a);
}

void slit_f(double *t, double *rV, double *F)
{
    *(F + 0) = *(rV + 2);
    *(F + 1) = *(rV + 3);
    *(F + 2) = 0;
    *(F + 3) = 0;
}

void slit_step(double *t, double *rV)
{
    *t += _dt;
    double X[4], k[4][4];
    for(int k_num = 0; k_num < 4; k_num++)
    {
        if (k_num == 0) for (int i = 0; i < 4; i++) X[i] = *(rV + i);
        if (k_num == 1) for (int i = 0; i < 4; i++) X[i] = *(rV + i) + k[0][i]/2.0;
        if (k_num == 2) for (int i = 0; i < 4; i++) X[i] = *(rV + i) + k[1][i]/2.0;
        if (k_num == 3) for (int i = 0; i < 4; i++) X[i] = *(rV + i) + k[2][i];
        double F[4];
        slit_f(t, rV, F);
        for (int i = 0; i < 4; i++)    k[k_num][i] = _dt*F[i];
    }
    for (int i = 0; i < 4; i++)
        *(rV + i) += (k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i])/6.0;
}

int slit_check(double *rV)
{
    if ((*(rV + 0) < 0) || (*(rV + 1) >= (_h + _dh)) || (*(rV + 1) <= (_h - _dh)))
        return 0;
    if (*(rV + 0) >= _l)
        return 2;
    return 1;
}

int slit_track(double *t, double *rV)
{
    slit_start(t, rV);
    while (1)
    {
        slit_step(t, rV);
        int dir = slit_check(rV);
        if (dir == 0)
            return 0;
        if (dir == 2)
            return 1;
    }
}

//=====================
//====    track    ====
//=====================

void track_f(double U, double *t, double *rV, double *F)
{
    *(F + 0) = *(rV + 2);
    *(F + 1) = *(rV + 3);
    *(F + 2) = 0;
    *(F + 3) = -(_qi/_mi)*(U/_H);
}

void track_step(double U, double *t, double *rV)
{
    *t += _dt;
    double X[4], k[4][4];
    for(int k_num = 0; k_num < 4; k_num++)
    {
        if (k_num == 0) for (int i = 0; i < 4; i++) X[i] = *(rV + i);
        if (k_num == 1) for (int i = 0; i < 4; i++) X[i] = *(rV + i) + k[0][i]/2.0;
        if (k_num == 2) for (int i = 0; i < 4; i++) X[i] = *(rV + i) + k[1][i]/2.0;
        if (k_num == 3) for (int i = 0; i < 4; i++) X[i] = *(rV + i) + k[2][i];
        double F[4];
        track_f(U, t, rV, F);
        for (int i = 0; i < 4; i++)    k[k_num][i] = _dt*F[i];
    }
    for (int i = 0; i < 4; i++)
        *(rV + i) += (k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i])/6.0;
}

int track_check(double *rV)
{
    if ((*(rV + 0) < 0) || (*(rV + 0) > (_L + _dL)) || (*(rV + 1) > _H))
        return 0;
    if (*(rV + 1) <= 0)
    {
        if ((*(rV + 0) >= (_L - _dL)) && (*(rV + 0) <= (_L + _dL)))
            return 2;
        return 0;
    }
    return 1;
}

int track_track(double U, double *t, double *rV)
{
    *t = 0;    
    while (1)
    {
        track_step(U, t, rV);
        int dir = track_check(rV);
        if (dir == 0)
            return 0;
        if (dir == 2)
            return 1;
    }
}

//=====================================
//========    Distributions    ========
//=====================================

double zero_one_dist()
{
    return 0.000001*(0 + rand() % 1000001);
}

double uniform(double start, double stop)
{
    return (start + (stop-start)*zero_one_dist());
}

double gauss(double E, double Var)
{
    double x, y, s;
    for(;;)
    {
        x = uniform(-1,1);
        y = uniform(-1,1);
        s = x*x + y*y;
        if (s>0 && s<=1)
            return (E + Var*x*sqrt((-2*log(s))/s));
    }
}
