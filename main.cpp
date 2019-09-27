#include <iostream>
#include <vector>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <map>
#include "FPCEC2010.h"
#include <algorithm>

using namespace std;


default_random_engine generator;
double* differential_evolution(int tamP, int dim , double CR, double F, int itMax, vector<double> lim1, vector<double> lim2, void (*cost_function)(double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh), int nh, int ng, double theta, double cp);
double phi_1(double* igualdad, double* desigualdad, int ng, int nh);
bool epsilon_comparacion(double f1, double f2, double phi1, double phi2, double epsilon);


int main()
{
    FILE* fp;
    fp=fopen("OTRO_eDE_C01.txt","wt");
    double mejor=(double)INT_MAX, peor=(double)INT_MIN, media=0, sd=0, v=0;
    int c=0;
    vector<double> soluciones;
    vector<int> violaciones;
    vector<double> phi;
    for(int it=0; it<1; it++)
    {
        cout<<it<<endl;
        int dim=10;
        vector<double> limL;
        vector<double> limU;
        for(int i=0; i<dim; i++)
        {
            limL.push_back(-1000);
            limU.push_back(1000);
        }
        double* resp=differential_evolution(200,dim,0.95,0.1,10000,limL,limU,C03,1,0,0.0001,1.0);
        soluciones.push_back(resp[dim]);
        phi.push_back(resp[dim+1]);
        violaciones.push_back(resp[dim+2]);
        printf("%e",soluciones[it]);
        cout<<" "<<phi[it]<<" "<<violaciones[it]<<endl;
    }
    for(int i=0; i<soluciones.size(); i++)
    {
        media+=soluciones[i];
        v+=phi[i];
        if(mejor>soluciones[i])
            mejor=soluciones[i];
        if(peor<soluciones[i])
            peor=soluciones[i];
        if(violaciones[i]==0)
            c++;
    }
    media=media/((double)soluciones.size());
    v=v/((double)phi.size());

    for(int i=0; i<soluciones.size(); i++)
        sd+=(soluciones[i]-media)*(soluciones[i]-media);
    sd=sqrt(sd/soluciones.size());
    sort(soluciones.begin(),soluciones.end());
    fprintf(fp,"%e %e %e %d %e %e %e",mejor,soluciones[(int)(soluciones.size()/2)],peor,c,v,media,sd);
    return 0;
}

double phi_1(double* igualdad, double* desigualdad, int ng, int nh)
{
    double sum=0;
    for(int j=0; j<nh; j++)
        if(fabs(igualdad[j])>=0.001)
            sum+=fabs(igualdad[j]);
    for(int j=0; j<ng; j++){
        double r = desigualdad[j];
        sum+=max(0.0, r);
    }
    return sum;
}

bool epsilon_comparacion(double f1, double f2, double phi1, double phi2, double epsilon)
{
    if(phi1<epsilon && phi2<epsilon)
    {
        if(f1<=f2)
            return true;
        else
            return false;
    }
    if(phi1==phi2)
    {
        if(f1<=f2)
            return true;
        else
            return false;
    }
    if(phi1<phi2)
        return true;
    else
        return false;
}

double* differential_evolution(int tamP, int dim , double CR, double F, int itMax, vector<double> lim1, vector<double> lim2, void (*cost_function)(double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh), int nh, int ng, double theta, double cp)
{
    //Generamos poblacion inicial
    int tamArchivo=2*tamP;
    double poblacion[tamArchivo][dim+3];
    uniform_real_distribution<double> distribution(0.0,1.0);
    double best[dim+3];

    FILE* ff;
    ff=fopen("eval.txt","w");
    FILE* ff1;
    ff1=fopen("eps.txt","w");

    //Generamos los valores de la Poblacion
    for(int i=0; i<tamArchivo; i++)
        for(int j=0; j<dim; j++)
            poblacion[i][j]=lim1[j]+(rand()%((int)(lim2[j]-lim1[j])))+ distribution(generator);

    int indice;
    double mejorF=999999999;
    double mejorPhi=999999999;

    //Evaluamos la aptitud de cada individuo y encontramos al mejor
    for(int i=0; i<tamArchivo; i++)
    {
        double aux[dim];
        for(int j=0; j<dim; j++)
            aux[j]=poblacion[i][j];

        double xx, ph[nh], pg[ng];
        cost_function(aux,&xx,pg,ph,dim,1,ng,nh);
        //Evaluamos la aptitud de cada individuo
        poblacion[i][dim]=xx;
        //Calculamos la funcion de violacion de restricciones
        poblacion[i][dim+1]=phi_1(ph,pg,ng,nh);
        int cont=0;
        for(int j=0; j <nh; j++)
            if(fabs(ph[j]) >= 0.01)
                cont++;
        for(int j=0; j <ng; j++)
            if(pg[j] >= 0.01)
                cont++;
        poblacion[i][dim+2]=cont;
    }
    double eps0, peor=0;
    for(int i=0; i<tamP; i++ )
        if(peor<poblacion[i][dim+1])
            peor=poblacion[i][dim+1];
    eps0=peor*theta;
    double eps=eps0;

    for(int i=0; i<tamArchivo; i++)
    {
        //Comparamos entre el nuevo individuo y el mejor hasta el momento
        bool comparacion= epsilon_comparacion(poblacion[i][dim],mejorF,poblacion[i][dim+1],mejorPhi,eps);
        if(comparacion==true)
        {
            mejorF=poblacion[i][dim];
            indice=i;
            mejorPhi=poblacion[i][dim+1];
        }
    }

    //Guardamos al mejor
    for(int i=0; i<=dim+2; i++){
        best[i]=poblacion[indice][i];
        poblacion[0][i]=best[i];
    }

    //Ciclo principal de la ED
    for(int i=0; i<itMax; i++)
    {
        cout<<"Mejor "<<i<<" -> ";
        printf("%e",best[dim]);
        fprintf(ff,"%e, ",best[dim]);
        fprintf(ff1,"%e, ",eps);
        cout<<" "<<best[dim+1]<<"-----"<<best[dim+2]<<endl;
        double new_poblacion[tamP][dim+3];
        for(int j=0; j<tamP; j++)
        {
            //Generamos 3 individuos nuevos y distintos entre si
            int r1, r2, r3;
            r1=rand()%tamP;
            while(r1==j)
                r1=rand()%tamP;
            r2=rand()%tamP;
            while(r2==j || r2==r1)
                r2=rand()%tamP;
            double bolado=distribution(generator);
            if(bolado<0.95)
                r3=tamP+rand()%tamP;
            else
                r3=rand()%tamP;
            while(r3==j || r3==r1 || r3==r2)
                r3=rand()%tamP;

            //Creamos la mutacion del nuevo individuo
            double ind[dim+3];
            for(int k=0; k<dim; k++)
                ind[k]=poblacion[r1][k] + F*(poblacion[r2][k]-poblacion[r3][k]);

            //Realizamos la cruza exponencial
            int kk, jj=rand()%dim;
            for(kk=0; kk<dim; kk++)
            {
                double rn=distribution(generator);
                if(!(rn<CR))
                    break;
                jj=(jj+1)%dim;
            }
            for(kk; kk<dim; kk++){
                ind[jj]=poblacion[j][jj];
                jj=(jj+1)%dim;
            }

            //Evaluamos el hijo creado
            double aux[dim];
            for(int k=0; k<dim; k++)
                aux[k]=ind[k];

            double xx, ph[nh], pg[ng];
            cost_function(aux,&xx,pg,ph,dim,1,ng,nh);
            //Evaluamos la aptitud de cada individuo
            ind[dim]=xx;
            //Calculamos la funcion de violacion de restricciones
            ind[dim+1]=phi_1(ph,pg,ng,nh);

            int cont=0;
            for(int t=0; t <nh; t++)
                if(fabs(ph[t]) >= 0.01)
                    cont++;
            for(int t=0; t<ng; t++)
                if(pg[t] >= 0.01)
                    cont++;
            ind[dim+2]=cont;

            //Torneo Binario usando la comparacion epsilon
            bool comparacion= epsilon_comparacion(ind[dim],poblacion[j][dim],ind[dim+1],poblacion[j][dim+1],eps);
            if(comparacion)
                for(int k=0; k<=dim+2; k++)
                    new_poblacion[j][k]=ind[k];
            else
                for(int k=0; k<=dim+2; k++)
                    new_poblacion[j][k]=poblacion[j][k];
        }

        //Copiamos la nueva poblacion
        for(int j=0; j<tamP; j++){
            for(int k=0; k<=dim+2; k++){
                poblacion[j][k]=new_poblacion[j][k];
            }
        }

        indice=-1;
        mejorF=9999999999;
        mejorPhi=9999999999;

        //Actualizamos el mejor
        for(int j=0; j<tamP; j++)
        {
            bool comparacion= epsilon_comparacion(poblacion[j][dim],mejorF,poblacion[j][dim+1],mejorPhi,eps);
            if(comparacion)
            {
                mejorF=poblacion[j][dim];
                indice=j;
                mejorPhi=poblacion[j][dim+1];
            }
        }
        bool comparacion= epsilon_comparacion(poblacion[indice][dim],best[dim],poblacion[indice][dim+1],best[dim+1],eps);
        if(comparacion)
            for(int j=0; j<=dim+2; j++)
                    best[j]=poblacion[indice][j];
        //Actualizar epsilon
        eps=eps0*pow((1.0 - ((double)i)/((double)itMax)),cp);
        cout<<"------------------------------->"<<eps<<"    "<<eps0<<endl;
    }

    return best;
}
