#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <conio.h>
#define clrscr() system("clear")
#define getch() system("read ")

#define eulernum 2.71828182846

typedef struct _fxpol
{
    double A;
    double B;
    double C;
    double D;
    double x;
} fxpol;

fxpol * derive_fx(fxpol * fxp, int veces)
{
    if(veces==0)
    	return fxp;
    else
    {
	fxp->A=fxp->A*fxp->B;
	fxp->B=fxp->B-1;
	fxp->C=fxp->C*(-1)*fxp->D;
	veces=veces-1;
	derive_fx(fxp,veces);
    }
}

/*void imp_exp(fxpol * fx,float A,float B,float C,float D)
{
    printf("Exp: %.0fx^%.0f+%.0fe^-%.0f\n",A,B,C,D);
    fx->A=1;fx->B=3;fx->C=1;fx->D=4;
}*/

fxpol * initfx(fxpol * fx)
{
    fx=malloc(sizeof(fxpol));
    fx->A=0;
    fx->B=0;
    fx->C=0;
    fx->D=0;
    return fx;
}

double eval_fx(double A,double B,double C,double D, double x)
{
    return A*pow(x,B)+C*pow(eulernum,(-1.0)*D*x);
}
int main()
{
    fxpol * fx1=initfx(fx1);
    fx1->A=1.0;fx1->B=3.0;fx1->C=1.0;fx1->D=4.0;fx1->x=2.0;
    printf("fx: %.30f\n",eval_fx(fx1->A,fx1->B,fx1->C,fx1->D,fx1->x));
    printf("df/dx: ");
    derive_fx(fx1,1);
    printf(" %.30f\n",eval_fx(fx1->A,fx1->B,fx1->C,fx1->D,fx1->x));
    return 0;
}

