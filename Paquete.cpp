/* LIBRERIAS */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define nqq 2
#define nqqq 3
#include <string.h>
#include <iostream>
using namespace std;
float e=2.71828;
/* DEFINICIONES */
#define PI 3.14159265

/* DECLARACIONES */
double *InicializarVector(int dim);
double **InicializarMatriz(int dim);

// Funcion necesarias para interpolar por lagrange
void Formula_Lagrange();

// Funcion necesaria para Interpolar por newton hacia adelante (Diferencias divididas)
void Newton_Adelante();

// Funcion necesaria para interpolar por newton hacia atras (Diferencias divididas)
void Newton_Atras();

void Spline();

int Integracion1();
int Integracion2();
double f1(double x);
double f2(double x);
double simpson13(float h, int n, double fx[50]);
double simpson38(float h, int n, int ultimos, double fx[50]);
double integracionNumerica(int f);

double jmf1, jmf2, jmf3;

	int inversa23(double mat[nqqq][nqqq],double FF,double FFF, double FFFF){

	int m=nqqq;
    int i, j, z, h, d;
	long double  a[m][m*2];


	long double matriz[m][m*2], dif, mult, f;
	for (i=0; i<m; i++) {
		for (j=0; j<m*2; j++) {
			matriz[i][j]=0;
		}
	}
	j=m;
	for (i=0; i<m; i++) {
		matriz[i][j]=1;
		j++;
	}

	for (i=0; i<m; i++) {
		for (j=0; j<m; j++) {
			 matriz[i][j]=mat[i][j];
		}
	}

	double cro=((matriz[0][0]*matriz[1][1]*matriz[2][2])+(matriz[0][1]*matriz[1][2]*matriz[2][0])+(matriz[0][2]*matriz[1][0]*matriz[2][1]));
	double cru=((matriz[2][2]*matriz[1][0]*matriz[0][1])+(matriz[2][1]*matriz[1][2]*matriz[0][0])+(matriz[2][0]*matriz[1][1]*matriz[0][2]));

	f=cro-cru;

	if(f>0 || f<0){

	for (z=0; z<m; z++) {
		dif=matriz[z][z];

		//Cambio de renglon
		if(dif==0){
			for(d=0;d<m*2;d++){
				a[0][d]=matriz[1][d];
				matriz[1][d]=matriz[0][d];
				matriz[0][d]=a[0][d];
			}

		}

		dif=matriz[z][z];

		for (h=0; h<m*2; h++) {
			matriz[z][h]=matriz[z][h]/dif;
		}
		for (i=0; i<m; i++) {
			if (i!=z){
				mult= -matriz[i][z];
				for (j=0; j<m*2; j++) {
					matriz[i][j]=matriz[z][j]*mult + matriz[i][j];
				}
			}
		}
	}

	}else{
		cout<<"\nNo se puede calcular la inversa de la matriz Jacobiana ya que el determinante es 0\n";
		return 1;
	}

	jmf1=matriz[0][3]*FF+matriz[0][4]*FFF+matriz[0][5]*FFFF;
	jmf2=matriz[1][3]*FF+matriz[1][4]*FFF+matriz[1][5]*FFFF;
	jmf3=matriz[2][3]*FF+matriz[2][4]*FFF+matriz[2][5]*FFFF;
	}

int inversa2(double mat[nqq][nqq],double FF,double FFF){

	int m=nqq;
    int i, j, z, h;
	double  a[m][m*2];
	double matriz[m][m*2], dif, mult, f;
	for (i=0; i<m; i++) {
		for (j=0; j<m*2; j++) {
			matriz[i][j]=0;
		}
	}
	j=m;
	for (i=0; i<m; i++) {
		matriz[i][j]=1;
		j++;
	}

	for (i=0; i<m; i++) {
		for (j=0; j<m; j++) {
			 matriz[i][j]=mat[i][j];
		}
	}

	double cro=(matriz[0][0]*matriz[1][1]);
	double cru=(matriz[1][0]*matriz[0][1]);

	f=cro-cru;

	if(f>0 || f<0){

	for (z=0; z<m; z++) {
		dif=matriz[z][z];

		//Cambio de renglon
		if(dif==0){
			for(int d=0;d<m*2;d++){
				a[0][d]=matriz[1][d];
				matriz[1][d]=matriz[0][d];
				matriz[0][d]=a[0][d];
			}

		}

		dif=matriz[z][z];

		for (h=0; h<m*2; h++) {
			matriz[z][h]=matriz[z][h]/dif;
		}
		for (i=0; i<m; i++) {
			if (i!=z){
				mult= -matriz[i][z];
				for (j=0; j<m*2; j++) {
					matriz[i][j]=matriz[z][j]*mult + matriz[i][j];
				}
			}
		}
	}

	}else{
		cout<<"\nNo se puede calcular la inversa de la matriz Jacobiana ya que el determinante es 0\n";
		return 1;
	}

	jmf1=matriz[0][2]*FF+matriz[0][3]*FFF;
	jmf2=matriz[1][2]*FF+matriz[1][3]*FFF;
	}
	
void newton(int op){

	int it;
	long double tol, Er, Ex, Ey, nx, ny, f1x, f1y, f2x, f2y, F1, F2, x, y, z, Ez, nz, f1z, f2z, f3x, f3y, f3z, F3;

	if(op>2){
		printf("\nDar el Punto Inicial (x, y, z)\n");
		printf("x: ");
		scanf("%f", &x);
		cout<<"y: ";
		cin>>y;
		cout<<"z: ";
		cin>>z;
	}else{
		cout<<"\nDar el Punto Inicial (x, y)\n";
		cout<<"x: ";
		cin>>x;
		cout<<"y: ";
		cin>>y;
	}

	cout<<"Dar el numero Maximo de Iteraciones:";
	cin>>it;
	cout<<"Dar el valor de Tolerancia:";
	cin>>tol;
	cout<<"\n";

	for(int freno=0;freno<it;freno++){

		cout<<"Iteracion "<<freno+1;
		switch(op){
			case 1:
				f1x=2*x+y;
				f1y=x;
				f2x=3*(y*y);
				f2y=6*x*y+1;
				F1=((x*x)+(x*y))-10;
				F2=(y+(3*x*(y*y)))-50;
				break;
			case 2:
				f1x=2*x;
				f1y=2*y;
				f2x=-exp(x);
				f2y=-2;
				F1=(pow(x, 2))+(pow(y, 2))-9;
				F2=-exp(x)-(2*y)-3;
				break;
			case 3:
				f1x=(4*x)-4;
				f1y=2*y;
				f1z=(6*z)+6;
				f2x=2*x;
				f2y=(2*y)-2;
				f2z=4*z;
				f3x=(6*x)-12;
				f3y=2*y;
				f3z=(-6*z);
				F1=((2*(pow(x,2)))-(4*x)+(pow(y,2))+(3*(pow(z,2)))+(6*z)+2);
				F2=((pow(x,2))+(pow(y,2))-(2*y)+(2*(pow(z,2)))-5);
				F3=((3*(pow(x,2)))-(12*x)+(pow(y,2))-(3*(pow(z,2)))+8);
				break;
			case 4:
				f1x=(2*x)-4;
				f1y=2*y;
				f1z=0;
				f2x=(2*x)-1;
				f2y=-12;
				f2z=0;
				f3x=(6*x)-12;
				f3y=2*y;
				f3z=(-6*z);
				F1=(pow(x,2)-(4*x)+(pow(y,2)));
				F2=((pow(x,2))-x-(12*y)+1);
				F3=((3*(pow(x,2)))-(12*x)+(pow(y,2))-(3*(pow(z,2)))+8);
				break;
		}

		if(op>2){
			double Jac[nqqq][nqqq];
			Jac[0][0]=f1x;
			Jac[0][1]=f1y;
			Jac[0][2]=f1z;
			Jac[1][0]=f2x;
			Jac[1][1]=f2y;
			Jac[1][2]=f2z;
			Jac[2][0]=f3x;
			Jac[2][1]=f3y;
			Jac[2][2]=f3z;

			if(inversa23(Jac,F1,F2,F3)==1){
			freno=it;
			}else{

			nx=x-jmf1;
			ny=y-jmf2;
			nz=z-jmf3;

			cout<<"\nPunto x= "<<x;
			cout<<"\nPunto y= "<<y;
			cout<<"\nPunto z= "<<z;

			cout<<"\n\nMatriz Jacobiana:\n";
			for (int i=0; i<nqqq; i++)
			{
				for (int j=0; j<nqqq; j++)
				{
					cout<<Jac[i][j];
					cout<<"   ";
				}
			cout<<"\n";
			}

			cout<<"\nf1(x,y,z)= "<<F1;
			cout<<"\nf2(x,y,z)= "<<F2;
			cout<<"\nf3(x,y,z)= "<<F3;

			cout<<"\n\nInversa de la matriz Jacobiana por f1= "<<jmf1;
			cout<<"\nInversa de la matriz Jacobiana por f2= "<<jmf2;
			cout<<"\nInversa de la matriz Jacobiana por f3= "<<jmf3;

			long double abs, absy, noy, nox, absz, noz;

			if(it>1){

			if((nx-x)<0){
				abs=(nx-x)*-1;
			}else{
				abs=nx-x;
			}
			if(nx<0){
				long double nox=nx*-1;
				Ex=abs/nox;
			}else{
				nox=nx;
				Ex=abs/nox;
			}

			if((ny-y)<0){

			long double absy=(ny-y)*-1;

			}else{
				absy=ny-y;
			}
			if(ny<0){
				long double noy=ny*-1;
				Ey=absy/noy;
			}else{
				noy=ny;
				Ey=absy/noy;
			}

			if((nz-z)<0){

			long double absz=(nz-z)*-1;

			}else{
				absz=nz-z;
			}
			if(nz<0){
				long double noz=nz*-1;
				Ez=absz/noz;
			}else{
				noz=nz;
				Ez=absz/noz;
			}

			if(Ex>Ey){
				Er=Ex;
			}else{
				Er=Ey;
			}if(Er>Ez){
				Er=Er;
			}else{
				Er=Ez;
			}

			cout<<"\n\nEr= "<<Er;

			if(Er<tol){
				cout<<"\nDetener";
				cout<<"\nEl metodo converge en el punto ("<<x<<", "<<y<<", "<<z<<") ";
				cout<<"en la iteracion: "<<freno+1<<" no es necesario hacer mas iteraciones\n\n\n";
				freno=it;
			}else{
				cout<<"\nSe deben hacer mas iteraciones ya que el metodo aun no converge\n\n\n";
			}

			}

			x=nx;
			y=ny;
			z=nz;
			}

		}else{
		double Jac[nqq][nqq];
		Jac[0][0]=f1x;
		Jac[0][1]=f1y;
		Jac[1][0]=f2x;
		Jac[1][1]=f2y;

		if(inversa2(Jac,F1,F2)==1){
			freno=it;
		}else{

		nx=x-jmf1;
		ny=y-jmf2;

		cout<<"\nPunto x= "<<x;
		cout<<"\nPunto y= "<<y;

		cout<<"\n\nMatriz Jacobiana:\n";
		for (int i=0; i<nqq; i++)
		{
			for (int j=0; j<nqq; j++)
			{
				cout<<Jac[i][j];
				cout<<"   ";
			}
		cout<<"\n";
		}

		cout<<"\nf1(x,y)= "<<F1;
		cout<<"\nf2(x,y)= "<<F2;

		cout<<"\n\nInversa de la matriz Jacobiana por f1= "<<jmf1;
		cout<<"\nInversa de la matriz Jacobiana por f2= "<<jmf2;

		long double abs, absy, noy, nox;

		if(it>1){


		if((nx-x)<0){
			abs=(nx-x)*-1;
		}else{
			abs=nx-x;
		}
		if(nx<0){
			long double nox=nx*-1;
			Ex=abs/nox;
		}else{
			nox=nx;
			Ex=abs/nox;
		}

		if((ny-y)<0){

		long double absy=(ny-y)*-1;

		}else{
			absy=ny-y;
		}
		if(ny<0){
			long double noy=ny*-1;
			Ey=absy/noy;
		}else{
			noy=ny;
			Ey=absy/noy;
		}

		if(Ex>Ey){
			Er=Ex;
		}else{
			Er=Ey;
		}

		cout<<"\n\nEr= "<<Er;

		if(Er<tol){
			cout<<"\nDetener";
			cout<<"\nEl metodo converge en el punto ("<<x<<", "<<y<<") ";
			cout<<"en la iteracion: "<<freno+1<<" no es necesario hacer mas iteraciones\n\n\n";
			freno=it;
		}else{
			cout<<"\nSe deben hacer mas iteraciones ya que el metodo aun no converge\n\n\n";
		}

		}

		x=nx;
		y=ny;

		}
		}

	}
}

void SE(int op){
	int resp;

	do{
		switch(op){
		case 1:
			cout<<"\nSE1\nf1(x, y)= x^2 + xy - 10 = 0\nf2(x, y)= y + 3xy^2 - 50 = 0";
			newton(op);
			break;
		case 2:
			cout<<"\nSE2\nf1(x, y)= x^2 + y^2 - 9 = 0\nf2(x, y)= -e^x - 2y - 3 = 0";
			newton(op);
			break;
		case 3:
			cout<<"\nSE3\nf1(x, y, z)= 2x^2 - 4x + y^2 + 3z^2 + 6z + 2 = 0\nf2(x, y, z)= x^2 + y^2 -2y + 2z^2 -5 = 0\nf3(x, y, z)= 3x^2 - 12x + y^2 -3z^2 + 8= 0";
			newton(op);
			break;
		case 4:
			cout<<"\nSE4\nf1(x, y, z)= x^2 + 4x + y^2 = 0\nf2(x, y, z)= x^2 -x -12y + 1 = 0\nf3(x, y, z)= 3x^2 -12x + y^2 -3z^2 + 8 = 0";
			newton(op);
			break;
		}
		cout<<"¿Desea probar con otros puntos? (1=Si, 2=No)";
		cin>>resp;
	}while(resp==1);
}

//DECLARACIONES GENERALES
double *InicializarVector(int dim)
{
    double *vector;
    vector = (double *)calloc(dim + 1, sizeof(double));
    if (vector == NULL)
    {
        printf("NO SE PUDO ASIGNAR MEMORIA\n\n");
        exit(0);
    }
    return vector;
}

double **InicializarMatriz(int dim)
{
    double **arreglo;
    arreglo = (double **)calloc(dim + 1, sizeof(double *));
    if (arreglo == NULL)
    {
        printf("\nNo se puede asignar memoria.\n");
        exit(1);
    }
    for(int m = 0; m <= dim + 1; m++)
    {
        arreglo[m] = (double *)calloc(dim + 1, sizeof(double));
        if (arreglo[m] == NULL)
        {
            printf("\nNo se puede asignar memoria.\n");
            exit(1);
        }
    }
    return arreglo;
}

int main(){
  	int n, opcion;
	int num, resp;
	int datos, op1;
  do{
        printf( "Bienvenido al paquete de M%ctodos Num%cricos II ", 130, 130);
        printf( "\nDesarrollado por: ");
        printf( "\n\tBuitron Arreola Juan Carlos");
        printf( "\n\tD%caz Garc%ca Flor", 161, 161);
        printf( "\n\tV%czquez P%crez X%cchitl", 160, 130, 162);
        printf( "\nA continuaci%cn las funciones del programa: \n", 162);
		    printf( "\n   1. Sistemas de Ecuaciones.");
        printf( "\n   2. Interpolaci%cn y suavizado de curvas.", 162);
        printf( "\n   3. Splines c%cbicos.", 163);
        printf( "\n   4. Diferenciaci%cn e Integraci%cn Num%crica.", 162, 162, 130);
        printf( "\n   5. Salir." );
        printf( "\n\n   Introduzca opci%cn (1-5): ", 162 );

        scanf( "%d", &opcion );

        /* Inicio del anidamiento */

        switch ( opcion )
        {
          case 1: do{
					system("cls");
					cout<<"METODO DE NEWTON-RAPHSON\n\n";
					cout<<"Programa hecho por Buitron Arreola Juan Carlos, Garcia Diaz Flor y Vazquez Perez Xochitl";
					cout<<"\nBuen dia, que Sistema de Ecuaciones no Lineales deseas resolver: ";
					cout<<"\n\nSE1\nf1(x, y)= x^2 + xy - 10 = 0\nf2(x, y)= y + 3xy^2 - 50 = 0";
					cout<<"\n\nSE2\nf1(x, y)= x^2 + y^2 - 9 = 0\nf2(x, y)= -e^x - 2y - 3 = 0";
					cout<<"\n\nSE3\nf1(x, y, z)= 2x^2 - 4x + y^2 + 3z^2 + 6z + 2 = 0\nf2(x, y, z)= x^2 + y^2 -2y + 2z^2 -5 = 0\nf3(x, y, z)= 3x^2 - 12x + y^2 -3z^2 + 8= 0";
					cout<<"\n\nSE4\nf1(x, y, z)= x^2 + 4x + y^2 = 0\nf2(x, y, z)= x^2 -x -12y + 1 = 0\nf3(x, y, z)= 3x^2 -12x + y^2 -3z^2 + 8 = 0";
					cout<<"\n\n5. Salir";
					cout<<"\n\nPresione el numero segun corresponda al SE que desea resolver [1-5]: ";
					cin>>num;
					system("cls");
					if(num<5){
					SE(num);
					cout<<"¿Desea resolver otro sistema? (1=Si, 2=No)";
					cin>>resp;
								}
				if(resp==2||num==5){
				cout<<"\n\nGracias por utilizar este programa buen dia";
				resp=2;
									  }
					}while(resp==1);
                    break;

            case 2: int Opcion;
                    system("cls");
                    printf("1.FORMULA DE LAGRANGE\n");
                    printf("2.FORMULA DE NEWTON (ADELANTE)\n");
                    printf("3.FORMULA DE NEWTON (ATRAS)\n");
                    printf("4.Salir\n");
                    printf("--------------------------------\n\n");
                    printf("Opcion: ");
                    scanf("%d", &Opcion);
                    switch (Opcion)
                    {
                    case 1:
                        system("cls");
                        Formula_Lagrange();
                        break;
                    case 2:
                        Newton_Adelante();
                        break;
                    case 3:
                        Newton_Atras();
                        break;
                    case 4:
                        exit(0);
                        break;
                    default:
                        printf("\n\nOpci%cn no encontrada\n\n\n", 162);
                    }
                    break;
            case 3:
                    Spline();
                    break;
            case 4:
                    Integracion1();
                    break;
            case 5:
                    exit(0);
                    break;
            default:
                    printf("\n\nOpci%cn no encontrada\n\n\n", 162);
        }
  }while(opcion!=1 || opcion!=2 || opcion!=3 || opcion!=4);
  
  return 0;
} 


//INTERPOLACION DE LAGRANGE
void Formula_Lagrange()
{
    int n;
    double *x, *y, x_i, suma = 0, prod = 0;
    printf("\n\nFORMULA DE LAGRANGE\n");
    printf("------------------------------------\n");
    printf("Ingrese numero de t%crminos: ", 130);
    scanf("%d", &n);
    printf("Ingrese el valor a interpolar: ");
    scanf("%lf", &x_i);
    x = InicializarVector(n);
    y = InicializarVector(n);
    printf("\nIngrese los datos de la tabla: \n");
    for (int i = 0; i < n; i++)
    {
        printf("x[%d] = ",i);
        scanf("%lf", &x[i]);
        printf("y[%d] = ",i);
        scanf("%lf",&y[i]);
    }
    printf("|  xi\t|yi     \n");
printf("--------------------\n");
    for (int i = 0; i < n; i++)
    {
        printf("|%.8lf\t|%.8lf\n", x[i], y[i]);
    }
    printf("\n");
    printf("Asi,\n");
    for (int i = 0; i < n; i++)
    {
        prod = 1.0;
        for (int j = 0; j < n; j++)
        {
            if (j != i)
            {
                prod *= ((x_i - x[j]) / (x[i] - x[j]));
                printf("[(%f - %f)/",x_i,x[j]);
                printf("(%f - %f)] ",x[i],x[j]);
            }
        }
        printf("%f ",y[i]);
        if(i+1<n){
             printf("+ \n");
        }
        suma += prod * y[i];
    }
    printf("\nFinalmente tenemos que...");
    printf("\n\n\Para x= %.8lf el valor de P(x) es: %.8lf\n\n",x_i, suma);
}

//INTERPOLACION DE NEWTON HACIA ADELANTE (DIFERENCIAS DIVIDIDAS)
void Newton_Adelante()
{
    int n = 0, i, k;                                       // Numero de suvintervalos
    double x_i = 0, *x, *y, **dif_div, suma = 0, prod = 0; // x_i es el punto a interpolar
    system("cls");
    printf("\n\t\t\tMETODO DE NEWTON\n");
    printf("\t\t (DIFERENCIAS DIVIDIDAS)\n");
    printf("\t...............................................\n");
    printf("\nIngrese el n%cmero de subintervalos: ", 163);
    scanf("%d", &n);
    printf("Ingrese el punto a interpolar: ");
    scanf("%lf", &x_i);
    x = InicializarVector(n);
    y = InicializarVector(n);
    for (i = 0; i <n; i++)
    {
       printf("x[%d] = ",i);
        scanf("%lf", &x[i]);
        printf("y[%d] = ",i);
        scanf("%lf",&y[i]);
    }
    system("cls");
    printf("\nMETODO DE NEWTON\n");
    printf("(DIFERENCIAS DIVIDIDAS)\n");
    printf("-----------------------------\n");
    printf("\nxi\t|yi\n");
    printf("-----------------------------\n");
    for (i = 0; i < n; i++)
    {
        printf("%lf|%.8lf\n", x[i], y[i]);
    }
    dif_div = InicializarMatriz(n);
    for (i = 0; i < n; i++)
    {
        if (i == 0)
        {
            for (k = 0; k < n; k++)
            {
                dif_div[k][i] = y[k];
            }
        }
        else
        {
            for (k = 0; k < n - i; k++)
            {
                dif_div[k][i] = (dif_div[k + 1][i - 1] - dif_div[k][i - 1]) / (x[k + i] - x[k]);
            }
        }
    }
    // Tabla de diferencias divididas
    printf("\n\n");
    for (i = 0; i < n; i++)
    {
        printf("\n");

        for (k = 0; k < n - i; k++)
        {
            printf("\t");
            printf("|%.8lf", dif_div[i][k]);
        }
    }
    suma = dif_div[0][0];
    for (i = 1; i < n; i++)
    {
        prod = 1.0;
        for (k = 0; k < i; k++)
        {
            prod *= (x_i - x[k]);
        }
        suma += dif_div[0][i] * prod;
    }
    printf("\n\nPara x= %.8lf, f(%.8lf)=%.8lf\n\n", x_i, x_i, suma);
}

//INTERPOLACION DE NEWTON HACIA ATRAS (DIFERENCIAS DIVIDIDAS)/
void Newton_Atras()
{
    int n = 0, i, k;                                       // Numero de suvintervalos
    double x_i = 0, *x, *y, **dif_div, suma = 0, prod = 0; // x_i es el punto a interpolar
    system("cls");
    printf("\n\t\t\tMETODO DE NEWTON\n");
    printf("\t\t (DIFERENCIAS DIVIDIDAS)\n");
    printf("\t...............................................\n");
    printf("\n\tIngrese el n%cmero de subintervalos: ", 163);
    scanf("%d", &n);
    printf("\tIngrese el punto 'x' a interpolar: ");
    scanf("%lf", &x_i);
    x = InicializarVector(n);
    y = InicializarVector(n);
    printf("\n\n\t\tx[i]\t\ty[i]\n");
    for (i = 0; i <= n; i++)
    {
        printf("\t\t");
        scanf("%lf\t\t%lf", &x[i], &y[i]);
    }
    system("cls");
    printf("\n\t\t\tMETODO DE NEWTON\n");
    printf("\t\t (DIFERENCIAS DIVIDIDAS)\n");
    printf("\t...............................................\n");
    printf("\n\t\tx[i]\t\t\ty[i]\n");
    printf("\t\t..................................\n");
    for (i = 0; i <= n; i++)
    {
        printf("\t\t");
        printf("%lf\t\t%.8lf\n", x[i], y[i]);
    }
    dif_div = InicializarMatriz(n);
    for (i = 0; i <= n; i++)
    {
        if (i == 0)
        {
            for (k = 0; k <= n; k++)
            {
                dif_div[k][i] = y[k];
            }
        }
        else
        {
            for (k = i; k <= n; k++)
            {
                dif_div[k][i] = (dif_div[k][i - 1] - dif_div[k - 1][i - 1]) / (x[k] - x[k - i]);
            }
        }
    }
    // Tabla de diferencias divididas
    printf("\n\n");
    for (i = 0; i <= n; i++)
    {
        printf("\n");
        for (k = 0; k <= i; k++)
        {
            printf("\t");
            printf("%.8lf", dif_div[i][k]);
        }
    }
    suma = dif_div[n][0];
    for (i = 1; i <= n; i++)
    {
        prod = 1.0;
        for (k = n; k > n - i; k--)
        {
            prod *= (x_i - x[k]);
        }
        suma += dif_div[n][i] * prod;
    }
    printf("\n\n\t\tPara x= %.8lf, f(%.8lf)=%.8lf\n\n", x_i, x_i, suma);
}
//Funciones Para Validar Enteros y Decimales
int ValidarEntero(){
	char numero[15];
	int validado, e;
	do{
		e=0;
		scanf("%s", &numero);
		for(int i=0; i<strlen(numero); i++){
		    if(!(isdigit(numero[i]))){
		        e++;
				if(i==0){
					if(numero[i]=='-'){
						e=e-1;
					}
				}
			}
		}
		if(e>0){
		        printf("\t\t***Error, ingresa un numero valido: ");
		}
	}while(e!=0);
	validado=atoi(numero);
	return validado;
}
float ValidarDecimal(){
	char numero[15];
	float validado;
	int e;
	do{
		e=0;
		scanf("%s", &numero);
		for(int i=0; i<strlen(numero); i++){
		    if(numero[i]!='.' && !(isdigit(numero[i])) ){
		        e++;
				if(numero[i]=='-'){
					e=e-1;
				}
			}
		}
		if(e>0){
		    printf("\t\t***Error, ingresa un n%cmero v%clido: ", 163, 160);
		}
	}while(e!=0);
	validado=atof(numero);
	return validado;
}
//SPLINE CUBICO
void Spline(){
	printf("SPLINE CUBICO METODO G_i(x)\n\n");
	int n;
	float x_inter;
	do{
		printf("\t-Cu%cntos Datos Vas A Usar, Contando x_0: ", 160);
		n=ValidarEntero();
		if(n<1){
			printf("\t***Error, solo se pueden ingresar datos positivos\n");
		}
	}while(n<1);
	int tam=n;
	float x[tam], y[tam], r[tam-1], dif[tam-1], h[tam-1];
	printf("\t-Ingresa los siguientes datos: \n");
	for(int i=0; i<tam; i++){
		printf("\t\tValor de x_%d: ", i);
		x[i]=ValidarDecimal();
		printf("\t\tValor de f(x_%d): ", i);
		y[i]=ValidarDecimal();
	}
	printf("\n\t- Ingresa El Valor x Que Quieres Interpolar: ");
	x_inter=ValidarDecimal();
	for(int i=0; i<(tam-1); i++){
		h[i]=x[i+1]-x[i];
		dif[i]=(y[i+1]-y[i])/h[i];
		r[i]=0;
	}
	float mat[tam-2][tam-1], s[tam], a[tam-1], b[tam-1], c[tam-1], d[tam-1];
	for(int j=0; j<(tam-2); j++){
		for(int k=0; k<(tam-1); k++){
			mat[j][k]=0;
		}
	}
	for(int i=0; i<(tam-2); i++){
		mat[i][tam-2]=(6*(dif[i+1]-dif[i]));
		if(i==0){
			mat[i][i]=2*(h[i]+h[i+1]);
			mat[i][i+1]=h[i+1];
		}else if(i==(tam-3)){
			mat[i][i-1]=h[i];
			mat[i][i]=2*(h[i]+h[i+1]);
		}else{
			mat[i][i-1]=h[i];
			mat[i][i]=2*(h[i+1]+h[i]);
			mat[i][i+1]=h[i+1];
		}
	}
	printf("\n\n|");
	for(int i=0; i<(tam-2); i++){
		printf("\t\t|");
		for(int j=0; j<(tam-2); j++){
			printf("   %f   ", mat[i][j]);
		}
		printf("| %f |\n", mat[i][tam-2]);
	}
	// Resolver Matriz
	float aux, pivoteo[tam-1];
	for(int i=0 ; i<(tam-2); i++){
		aux=mat[i][i];
		for(int j=0 ; j<(tam-1); j++){
			if(mat[i]!=0){
				mat[i][j]=mat[i][j]/aux;
				pivoteo[j]=mat[i][j] ;
			}
		}
		if(i<(tam-3)){
			aux=mat[i+1][i];
			for(int j=i; j<(tam-1); j++){
				mat[i+1][j]=mat[i+1][j]-(pivoteo[j]*aux);
			}
		}
	}
	for(int i=(tam-3); i>0; i--){
		aux=mat[i-1][i];
		for(int j=0 ; j<(tam-1); j++){
			if(mat[i]!=0){
				pivoteo[j]=mat[i][j] ;
			}
		}
		for(int j=0; j<(tam-1); j++){
				mat[i-1][j]=mat[i-1][j]-(pivoteo[j]*aux);
		}
	}
	printf("\n\n");
	for(int i=0; i<(tam-2); i++){
		printf("\t\t|");
		for(int j=0; j<(tam-2); j++){
			printf("   %f   ", mat[i][j]);
		}
		printf("| %f |\n", mat[i][tam-2]);
	}
	//
	printf("\n\n\t-Las ecuaciones resultantes son: \n");
	s[0]=0, s[tam-1]=0;
	for(int i=0; i<tam; i++){
		if(i==0 || i==(tam-1)){
			s[i]=0.000000;
		}else{
			s[i]=mat[i-1][tam-2];
		}
	}
	int j;
	for(int i=0; i<(tam-1); i++){
		j=i+1;
		a[i]=(s[j]-s[i])/(6*h[i]);
		b[i]=s[i]/2;
		c[i]=dif[i]-(h[i]*((s[j]+(2*s[i]))/6));
		d[i]=y[i];
	}
	for(int i=0; i<(tam-1); i++){
			printf("\t\tg_%d(x) = (%f)(x-%f)^3 + (%f)(x-%f)^2 + (%f)(x-%f) + %f \n", i, a[i], x[i], b[i], x[i], c[i], x[i], d[i]);
	}
	float xs;
	printf("\n\n\tEl valores de los polinomios en el punto quieres interpolar son: \n");
	for(int i=0; i<(tam-1); i++){
		xs=x_inter-x[i];
		r[i]=r[i]+(a[i]*(pow(xs,3)));
		r[i]=r[i]+(b[i]*(pow(xs,2)));
		r[i]=r[i]+(c[i]*xs);
		r[i]=r[i]+d[i];
		printf("\t\tg_%d(%x_inter) = %f \n", i, x_inter, r[i]);
	}
	system("pause & cls");
	return;
}

//INTEGRACION
double f1(double x){
	double f1;	
		f1=(pow(x,4)*sqrt(3+2*pow(x,2)))/3;
	return f1;
}

double f2(double x){
	double f2;	
		f2=pow(x,5)/pow(pow(x,2)+4,0.2);
	return f2;
}

double simpson13(float h, int n, double fx[50]){
	int i;
	
	double simp13=(fx[0]+fx[n])*h/3;	
	
	for(i=1; i<n; i++){
		if(i%2==0){
			simp13=simp13+(2*fx[i]*h/3);
		}
		if(i%2!=0){
			simp13=simp13+(4*fx[i]*h/3);
		}
	}
	
	return simp13;
}

double simpson38(float h, int n, int ultimos, double fx[50]){
	int i;
	double simp38=(3*fx[n-ultimos]*h/8)+(3*fx[n]*h/8);
	
	
	for(i=n-ultimos+1; i<n; i++){
		simp38=simp38+(3*h*3*fx[i]/8);
	}
	
	return simp38;
}

double integracionNumerica(int f){
	int n=0, i=0;
	bool ok=true;
	double h, xini, xfin;
	double fx[50]={0}, solucion=0;
	
	do{
		ok=true;
		cout<<"\n\n\t\tIntervalo de integracion";
		cout<<"\n\n\tPunto de inicio: ";
		cin>>xini;
		cout<<"\n\n\tPunto final: ";
		cin>>xfin;
		
		if(xini>xfin){
			cout<<"\n\n\t\tEscriba el intervalo con el punto de inicio menor al punto final";
			ok=false;
		}
	}while(ok==false);
	
	do{
		ok=true;
		cout<<"\n\n\tNumero de subintervalos deseados ";
		cin>>n;
		if(n>50){
			cout<<"\n\n\tEl maximo de subintervalos permitidos es 50";
			ok=false;
		}
	}while(ok=false);
	
	h=(xfin-xini)/n;
	
	for(i=0; i<n+1; i++){
		if(f==1){
			fx[i]=f1(xini+(h*i));
		}
		if(f==2){
			fx[i]=f2(xini+(h*i));
		}
	}
	
	if(n%2==0){
		solucion=simpson13(h, n, fx);

	}else{
		solucion=simpson13(h, n-3, fx);
		solucion=solucion+simpson38(h, n, 3, fx);
		cout<<"\n\t"<<solucion;
	}
	
	return solucion;
}

int Integracion2(){
	int opcion;
  
  cout<<"\n\nSeleccione la funcion a integrar: \n\n", 162;
	
	cout<<"\t1. f(x)=(x^4)((3+2x^2)^1/2)/3)\n";
	cout<<"\t2. f(x)=(x^5)/((x^2+4)^1/5)\n";
	
	cout<<"\n\n\t3. Salir";

	do{
		cout<<"\n\n\t\tOpcion: ";
		cin>> opcion;
	}while(opcion<1 || opcion>3);
	
	return opcion;	
}
int Integracion1(){
	int salir=0, opcion;
	double solucion;
	
	do{  
	solucion=0;
	system("cls");
	
		opcion=Integracion2();
			
		switch(opcion){
			case 1:
				solucion = integracionNumerica(opcion);	
				cout<<"\n\tLa solucion es: "<<solucion;
				cout<<"\n\n\n";	
				system("pause");						
				break;
			case 2:
				solucion = integracionNumerica(opcion);	
				cout<<"\n\tLa solucion es: "<<solucion;	
				cout<<"\n\n\n";	
				system("pause");						
				break;
			case 3:
				return salir=1;
			}
	}while(salir!=1); 
}
