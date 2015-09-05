#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <conio.h>
#define clrscr() system("clear")
#define getch() system("read ")

#define eulernum 2.71828182846

typedef struct _fxpol
{
    float A;
    float B;
    float C;
    float D;
    float x;
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

void imp_exp(fxpol * fx,float A,float B,float C,float D)
{
    printf("Exp: %.0fx^%.0f+%.0fe^-%.0f\n",A,B,C,D);
}

fxpol * initfx(fxpol * fx)
{
    fx=malloc(sizeof(fxpol));
    fx->A=0;
    fx->B=0;
    fx->C=0;
    fx->D=0;
    return fx;
}

int main()
{
    fxpol * fx1=initfx(fx1);
    fx1->A=1;fx1->B=3;fx1->C=1;fx1->D=4;
    imp_exp(fx1,fx1->A,fx1->B,fx1->C,fx1->D);
    derive_fx(fx1,1);
    printf("df/dx: ");
    imp_exp(fx1,fx1->A,fx1->B,fx1->C,fx1->D);
    return 0;
}
