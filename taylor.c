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
double fxx=0;
double eabs=0;
double erel=0;
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

void initz()
{
    fact=1.0;
    tayl=0;
    fxx=0;
    eabs=0;
    erel=0;    
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
    printf("Exp: %.0fx^%.0f+%.0fe^-%.0fx\n",fx->A,fx->B,fx->C,fx->D);
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

double ABSn(double err)
{
    if(err<0)
        return -err;
    else
        return err;
}

double er_abs(double tayl)
{
    double err=ABSn(fxx-tayl);
    return err;
}

double er_rel(double tayl)
{
    double err=2*(er_abs(tayl)/(ABSn(tayl)+ABSn(fxx)));
    return err;
}


double eval_fx(double A,double B,double C,double D, double x)
{
    if (B>=0)
    	return A*pow(x,B)+C*pow(eulernum,(-1.0)*D*x);
    else
    	return A+C*pow(eulernum,(-1.0)*D*x);
}


double taylors(fxpol * fx,double ic)
{
     double err_ABS=er_abs(tayl);
     double err_REL=er_rel(tayl);
     printf("EABS: %.30f\n",err_ABS);
     printf("EREL: %.30f\n",err_REL);
     if(err_ABS<=epsilon)
     {
	fact=1;
	tayl+=((eval_fx(fx->A,fx->B,fx->C,fx->D,fx->x0))*pow((fx->x-fx->x0),ic))/factorial(ic);
        printf("IC: %.0f;\tTaylor:\t%.30f\n",ic,tayl);
	ic=0;
	eabs=err_ABS;
	erel=err_REL;
	return tayl;
     }
     else
     {
	fact=1;
	tayl+=((eval_fx(fx->A,fx->B,fx->C,fx->D,fx->x0))*pow((fx->x-fx->x0),ic))/factorial(ic);
        printf("IC: %.0f;\tTaylor:\t%.30f\n",ic,tayl);
	ic++;
	derive_fx(fx,1);
	taylors(fx,ic);
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
//TODO: La derivada al evaluarla con 0 causa estragos, no se está haciendo la derivada real, se está simulando.
//Parcialmente reparada la derivada.

void ntayl()
{
    clrscr();
    initz();
    fxpol * fx1=initfx(fx1);
    //fx1->A=1;fx1->B=3;fx1->C=1;fx1->D=4;fx1->x=2;fx1->x0=1;
    printf("Ingrese los coeficientes para la función de tipo: Ax^B+Ce^(-Dx)\n");
    input(fx1);
    fxx=eval_fx(fx1->A,fx1->B,fx1->C,fx1->D,fx1->x);
    printf("Resultado de la serie de Taylor: %.5f\n",taylors(fx1,0));
    printf("Error Absoluto: %.5f\n",eabs);
    printf("Error Relativo Porcentual: %.5f%\n",erel*100);
    printf("Valor Exacto de la función: %.5f\n",fxx);
    getch();
}
void presentacion()
{
	clrscr();
	printf("\n                      Universidad Tecnol¢gica de Panama \n");
	printf("            Facultad de Ingenier¡a en Sistemas Computacionales\n");
	printf("           Departamento de Computaci¢n y Simulaci¢n de Sistemas\n\n");
	printf("                         Métodos Numéricos para Ingenieros\n\n");
	printf("                               Proyecto nø1\n");
	printf("                        Ingeniera Jacqueline De Ching\n\n");
	printf("                                   Grupo: \n");
	printf("                          Espinosa, Angel  8-905-1352    \n"); 
    printf("                             Miranda, Yaneys  8-879-525\n\n");
	printf("                                   1IL-121\n");
	getch();
}

void menup()
{
    int op, cont=1, lmp=1;
    char * opp=malloc(sizeof(char));
    do
    {
	if (lmp)
            clrscr();
	printf("\t\t\tMENU PRINCIPAL\n");
	printf("1. Presentacion\n");
	printf("2. Cálculo de la Serie de Taylor\n");
	printf("3. Salir del Programa\n");
	printf("Ingrese una opcion: ");
	scanf("%s",opp);
	op=atoi(opp);
	switch (op)
	{
	    case 1: presentacion();		    
		    lmp=1;
		    break;
	    case 2: ntayl(); break;
	    case 3: cont=0; break;
	    default: printf("Escriba una opcion correcta.\n"); break;
	}
    } while( cont );
}



int main()
{
    menup();
    return 0;
}
