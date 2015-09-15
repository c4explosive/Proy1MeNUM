#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <conio.h> //Descomenta esta, para usarlo en turbo C.
////////////////////////////////////////////
#define clrscr() system("clear") //Comenta estas dos líneas
#define getch() system("read ")
////////////////////////////////////////////

#define eulernum 2.71828182845904 //Constante de euler
#define epsilon 0.0001 //valor donde el error relativo tiene que llegar

double fact=1.0; //Variable de de factorial
double tayl=0; //Se debe inicializar para evitar sorpresas
typedef struct _fxpol //estructura del polinomio Ax^B+Ce^-Dx; donde A,B,C y D son coeficientes
{
    double A;
    double B;
    double C;
    double D;
    double x; //esta varible, n y x0 no son parte del polinomio pero si para evaluar la función al momento de leerla
    double x0;
    int n;
} fxpol; //creamos el tipo de dato que representa el polinomio

fxpol * derive_fx(fxpol * fxp, int veces) //Aquí es donde se deriva, lo hice recursivo, pero en realidad en taylor siempre se pone una vez, ya que esto se va acumulando, si le pongo que haga 3 veces es como si le dijera que en taylor derivara tres veces en un sólo término.
{
    if(veces==0) //Proceso recursivo, termina cuando las veces a derivar sean cero y retorna los nuevos valores de coeficientes.
    	return fxp;
    else
    {
	fxp->A=fxp->A*fxp->B; //Aquí básicamente no se deriva, realmente si no que cambia los valores de los coefientes, lo que causa el efecto de derivar.
	fxp->B=fxp->B-1;
	fxp->C=fxp->C*(-1)*fxp->D;
	veces=veces-1;
	derive_fx(fxp,veces);
    }
}
double factorial(double x) //Bueno el factorial por definición pero de forma recursiva n*fact(n-1)
{
     if(x==0)
	return fact*1;
     else
	fact=x*factorial(x-1);
}
void imp_pol(fxpol * fx) //Función de relleno, xD; era para ver como iba cambiando el polinomio a medida que deribava.
{
    printf("Exp: %.0fx^%.0f+%.0fe^-%.0fx\n",fx->A,fx->B,fx->C,fx->D);
    //fx->A=1;fx->B=3;fx->C=1;fx->D=4;
}

fxpol * initfx(fxpol * fx) //Esta función inicializa los coeficientes del polinomio para que no tengan basura.
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

void input(fxpol * fx)//Admito que esto es un enredo -_- xD... Pero es una manera de leer los datos por medio de punteros para no repetir la función...
{
    char Vars[][3]={"A","B","C","D","x","x0"}; //Esta array sirve para mostrar el nombre del elemento al usuario.
    int i; 
    double * val=NULL; //Este puntero es el que iremos aumentando para que cambie a la dirección de memoria a la que se quiere meter el dato
    val=&(fx->A); //Aquí el puntero toma la dirección base del primer elemento 
    for(i=0;i<6;i++) //Y en este ciclo va aumentando el valor de la dirección de modo que vaya ubicando los datos como si se tratase de un arreglo
    {
    	printf("Ingese el valor de %s: ",Vars[i]);
	scanf("%le",val);
	val++; //Aumenta la dirección de memoria
    }

}

double eval_fx(double A,double B,double C,double D, double x) //Esta función es donde se evalua el polinomio con todos los coeficientes y valor de x; retorna el resultado de forma flotante
{
    if (B>=0) //Esto es para cuando la derivada o el exponente de x sea 0, si esto sucede la expresión si se transforma y si lo dejo evaluando normal va causar una 1/0 y un agujero negro en Panamá...
    	return A*pow(x,B)+C*pow(eulernum,(-1.0)*D*x);
    else
    	return A+C*pow(eulernum,(-1.0)*D*x);
}


double taylors(fxpol * fx,int terms,double ic) //Bueno, aquí esta Taylor, tambien es recursivo y funciona así: recive el polinomio de la estructura, la cantidad de terminos a evaluar y la var ic es la que dice cuantas veces a iterado. Esta iteración es progresiva y tiene que tener la var terms que le dice hasta donde ic debe llegar y retornar el resultado.
{
     if(ic==terms) //Aquí es donde termina, se pone el fact en 1, se suma a la var tayl toda la formula de taylor y luego se retorna, aquí ya no se deriva.
     {
	fact=1;
	tayl+=((eval_fx(fx->A,fx->B,fx->C,fx->D,fx->x0))*pow((fx->x-fx->x0),ic))/factorial(ic);
        printf("IC: %.0f;\tTaylor:\t%.30f\n",ic,terms,tayl);
	ic=0;
	return tayl;
     }
     else  //bueno, no es experimentado pero creo queaquí uno de esos else está de más, ya que los dos derivan y hacen los mismo, aquí se define el factorial como 1, se suma a taylor un termino de su serie y se evalua ya así hasta que llegue al caso base.
     {
	fact=1;
	tayl+=((eval_fx(fx->A,fx->B,fx->C,fx->D,fx->x0))*pow((fx->x-fx->x0),ic))/factorial(ic);
        printf("IC: %.0f;\tTaylor:\t%.30f\n",ic,terms,tayl);
	ic++;
	derive_fx(fx,1);
	taylors(fx,terms,ic);
     }
}

void print_terms(fxpol * fx) //Otra función de relleno... Coma para visualizar el estado de cada coeficiente, funciona de la misma forma que al leer los datos.
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
int main() //Y aquí comienza el programa, se crea el polinomio y se inicializa con initfx() para que no tenga basura.
{
    fxpol * fx1=initfx(fx1);
    int n=500; //Esto se debe leer, es la cantidad de términos que se quiere de la serie
    fx1->A=1;fx1->B=3;fx1->C=1;fx1->D=4;fx1->x=2;fx1->x0=1; //Esta linea debe ser reemplazada por el input de abajo, era sólo de prueba
    //input(fx1);
    printf("fx: %.30f\n",eval_fx(fx1->A,fx1->B,fx1->C,fx1->D,fx1->x)); //Imprime el valor de la función, evaluado directaaamente.
    printf("Taylors Series: %.30f\n",taylors(fx1,n-1,0)); //Imprime la serie de taylor.
    /*derive_fx(fx1,4);
    print_terms(fx1);
    printf("fx: %.30f\n",eval_fx(fx1->A,fx1->B,fx1->C,fx1->D,fx1->x0));*/
    imp_pol(fx1);//Otro artilugio más para verificar el estado del polinomio
    return 0; //Como no uso turbo C, se debe retornar un estado al SO con esta línea.
}


