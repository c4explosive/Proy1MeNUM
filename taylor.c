#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <conio.h>
#define clrscr() system("clear")
#define getch() system("read ")

#define eulernum 2.71828182845904
#define epsilon 0.0001

double fact=1.0;
double tayl=0;
typedef struct _fxpol
{
    double A;
    double B;
    double C;
    double D;
    double x;
    double x0;
    int n;
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
double factorial(double x)
{
     if(x==0)
	return fact*1;
     else
	fact=x*factorial(x-1);
}
void imp_pol(fxpol * fx)
{
    printf("Exp: %.0fx^%.0f+%.0fe^-%.0f\n",fx->A,fx->B,fx->C,fx->D);
    //fx->A=1;fx->B=3;fx->C=1;fx->D=4;
}

fxpol * initfx(fxpol * fx)
{
    fx=malloc(sizeof(fxpol));
    fx->A=0;
    fx->B=0;
    fx->C=0;
    fx->D=0;
    fx->x=0;
    fx->x0=0;
    fx->n=0;
    return fx;
}

void input(fxpol * fx)
{
    char Vars[][3]={"A","B","C","D","x","x0"};
    int i;
    double * val=NULL;
    val=&(fx->A);
    for(i=0;i<6;i++)
    {
    	printf("Ingese el valor de %s: ",Vars[i]);
	scanf("%le",val);
	val++;
    }

}
/*double ipow(double b, double n)
{
    printf("N: %f\n",n);
    if (n<0)
	return (1/pow(b,-n));
    else
	return pow(b,n);
}*/

double eval_fx(double A,double B,double C,double D, double x)
{
    return A*pow(x,B)+C*pow(eulernum,(-1.0)*D*x);
}


double taylors(fxpol * fx,int terms,double ic)
{
     if(ic==terms)
     {
	fact=1;
	tayl+=((eval_fx(fx->A,fx->B,fx->C,fx->D,fx->x0))*pow((fx->x-fx->x0),ic))/factorial(ic);
        printf("IC: %.0f;\tTaylor:\t%.30f\n",ic,terms,tayl);
	ic=0;
	return tayl;
     }
     else if (ic==0)
     {
	fact=1;
	tayl+=((eval_fx(fx->A,fx->B,fx->C,fx->D,fx->x0))*pow((fx->x-fx->x0),ic))/factorial(ic);
        printf("IC: %.0f;\tTaylor:\t%.30f\n",ic,terms,tayl);
	ic++;
	derive_fx(fx,1);
	taylors(fx,terms,ic);
     }
     else
     {
	fact=1;
	tayl+=((eval_fx(fx->A,fx->B,fx->C,fx->D,fx->x0))*pow((fx->x-fx->x0),ic))/factorial(ic);
        printf("IC: %.0f;\tTaylor:\t%.30f\n",ic,terms,tayl);
	ic++;
	derive_fx(fx,1);
	taylors(fx,terms,ic);
     }
}

void print_terms(fxpol * fx)
{
    double * var=NULL;
    char Vars[][3]={"A","B","C","D","x","x0"};
    int i;
    var=&(fx->A);
    for(i=0;i<6;i++)
    {
	printf("%s: %.2f\n",Vars[i],*var);
	var++;
    }
}
//TODO: La derivada al evaluarla con 0 causa estraagos, no se está haciendo la derivada real, se está simulando.
int main()
{
    fxpol * fx1=initfx(fx1);
    int n=7;
    fx1->A=1;fx1->B=3;fx1->C=1;fx1->D=4;fx1->x=2;fx1->x0=1;
    //input(fx1);
    printf("fx: %.30f\n",eval_fx(fx1->A,fx1->B,fx1->C,fx1->D,fx1->x));
    printf("Taylors Series: %.30f\n",taylors(fx1,n-1,0));
    /*derive_fx(fx1,4);
    print_terms(fx1);
    printf("fx: %.30f\n",eval_fx(fx1->A,fx1->B,fx1->C,fx1->D,fx1->x0));*/
    return 0;
}


