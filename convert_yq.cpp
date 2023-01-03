//
//  NStoHIC
//
//  Created by Nanxi Yao on 2/26/22.
//


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <math.h>

using namespace std;

const double nsat =0.161;
const double hc = 197.3;
const double mN=938.0;
const double hbar=6.582119569*pow(10.0,-22.0);
const double c=2.998*pow(10.0,23.0);
const double mmu=105.7/c/c;
const double me=0.51099895/c/c;

const double PI  =3.141592653589793238463;

double getFermiEnergy(double mass, double ni)
{
    double kF=pow(ni/2.0*6*PI*PI,1.0/3.0)*hbar;
    double epsilon=pow(mass,4.0)*pow(c,5.0)/PI/PI/hbar/hbar/hbar;
    double x=kF/mass/c;
    double eden=epsilon/8.0*((2*pow(x,3.0)+x)*pow((1+x*x),1.0/2.0)-asinh(x));
    return eden;
    
}
double defYQ(double nB, double E, double L, double K)
{
    double x=(nB-nsat)/3.0/nsat;
    double power=(E+L*x+K*x*x/2.0)/hc;
    double yq=64.0/3.0/PI/PI/nsat/(3.0*x+1.0)*power*power*power;
    double powernsat = E/hc;
 //   double yqnsat =64.0/3.0/PI/PI/nsat*powernsat*powernsat*powernsat;
    if(yq>1)  return 1.0;
    if(yq<0)  return 0.0;
    else
        return yq;
}
/*Define a smooth yQ expansion which converges to zero*/
double getSmoothfunc(double yqval, double switch_pt, double x)
{
    return yqval/switch_pt*x;
}


double YQ_cvgt(double switch_nB, double E, double L, double K, double nB)
{
    double smoothness = 0.025;
    double g=defYQ(nB, E, L, K);
    double f=getSmoothfunc(defYQ(switch_nB,E, L, K),switch_nB,nB);
    double s=0.5+0.5*tanh((nB-switch_nB)/smoothness);
    double yQ = s*g+(1-s)*f;
    if (yQ>0.3) return 0.3;
    else
        return yQ;
}

// This conversion only considers expanding from NS to symmetric matter
double getEHIC(double eNS,double rho,double e0,double L, double K, double J, double yQ)
{
    double delta=0;
    return eNS-(e0+L*(rho/nsat-1)/3+K/18*(rho/nsat-1)*(rho/nsat-1)+J/162*((rho/nsat-1)*(rho/nsat-1)*(rho/nsat-1)))*(1-2*yQ)*(1-2*yQ)*rho;
}

//expanding from NS to symmetric matter and then to HIC (y_HIC varies)
double getEHIC2(double eNS,double rho,double e0,double L, double K, double J, double Y_HIC, double switch_pt)
{
  //  double delta=0;
    double yQ=YQ_cvgt(switch_pt, e0, L, K, rho);
  //  delta = 1-2*Y_HIC;
    return eNS-(e0+L*(rho/nsat-1)/3+K/18*(rho/nsat-1)*(rho/nsat-1)+J/162*((rho/nsat-1)*(rho/nsat-1)*(rho/nsat-1)))*4*((Y_HIC-yQ)+(yQ*yQ-Y_HIC*Y_HIC))*rho;
 //   return eNS-(e0+L*(rho/nsat-1)/3+K/18*(rho/nsat-1)*(rho/nsat-1)+J/162*((rho/nsat-1)*(rho/nsat-1)*(rho/nsat-1)))*delta*delta;
}

//checks if converted cs2 is causal for nB<rholimit, number of acausal points less than numlimit
bool isCausal(double *cs2, double rhoB[], double rholimit, double rholimit2, double numlimit, int limit)
{
    int count=0;
    for(int i=0;i<limit;i++)
    {
        if (rhoB[i]/nsat>rholimit && rhoB[i]/nsat<rholimit2)
        {
            
            if(cs2[i]<-0.03 || cs2[i]>1)
            {
                count++;
              //  cout<<cs2[i]<<endl;
            }
         //   cout<<cs2[i]<<" "<<count<<endl;
        }
    }
    if(count > numlimit)
    { //cout<<"num "<<count<<endl;
        return false;
    }
    return true;
    
}


double getnumericalderiv(double x[], double y[], int i)
{
    double deriv = (y[i+1]-y[i-1])/(x[i+1]-x[i-1]);
    return deriv;
}


double** convert(double *&rhoB, double *&energy, double e0, double L, double K, double J, int limit, double YHIC, double switch_pt)
{
    double enerHIC[limit], eden[limit], pressure[limit],cs2[limit],eoverrho[limit],rhoB2[limit];
    double *ne =new double[limit];
    double *energysub=new double[limit];
    for(int i=0;i<limit;i++)
    {
        double yQ=YQ_cvgt(switch_pt, e0, L, K, rhoB[i]);
        ne[i]=rhoB[i]*yQ;
        if(ne[i]<0) ne[i]=-ne[i];
        double neval=ne[i];
        energysub[i]=energy[i]-getFermiEnergy(me,neval);
      
    }
    
    for(int i=0;i<limit;i++)
    {
        eden[i]=getEHIC2(energysub[i],rhoB[i],e0,L,K,J, YHIC,switch_pt);
       // eden[i]=getEHIC(energy[i],rhoB[i],e0,L,K,J,YQ[i]);
        eoverrho[i]=eden[i]/rhoB[i];
        rhoB2[i]=rhoB[i];
    }
    for(int i=0;i<limit;i++)
    {
        pressure[i]=getnumericalderiv(rhoB2,eoverrho,i)*rhoB2[i]*rhoB2[i];
    }

    for(int i=0;i<limit;i++)
    {
        cs2[i]=getnumericalderiv(eden,pressure,i);
    }
    double** rhoB_cs2 = 0;
    rhoB_cs2 = new double*[2];
    rhoB_cs2[0] = new double[limit];
    rhoB_cs2[1] = new double[limit];
        
    for (int w = 0; w < limit; w++)
        {
        rhoB_cs2[0][w]=rhoB2[w];
     //   rhoB_cs2[1][w]=energysub[w];
    //  rhoB_cs2[1][w]=fermi[w];
      rhoB_cs2[1][w]=cs2[w];
   //   cout<<w<<" "<<rhoB2[w]<<endl;
    //  rhoB_cs2[1][w]=pressure[w];
    //  rhoB_cs2[1][w]=eden[w];
    //  rhoB_cs2[1][w]=eoverrho[w];
                }
  /*  for(int i=0;i<limit;i++)
    {
    //    cout<<i<<" "<<cs2<<endl;
    }*/
    return rhoB_cs2;
  //  return cs2;
}


int main() {
 
    FILE *finput = fopen("eos3_2.txt","r");
    int limit=6075;             //**Open the input file.
    double *energy = new double[limit];
    double *muB = new double[limit];
    double *rhoB = new double[limit];
    double *pressure =new double[limit];
    double drop2;

    for(int i=0;i<limit;i++)
    {
      fscanf(finput,"%lf\t%lf\t%lf\t%lf\t%lf\n",&rhoB[i],&pressure[i],&energy[i],&drop2,&muB[i]);
    }
    fclose(finput);
    /*get original EOS info*/
 /*   ofstream myfile2("original_eos4_cs2.txt");
   //  ofstream myfile2("sly4_allyq.txt");
  //
    double oldcs2[limit];
    double allyQ[limit];
    for(int i=0;i<limit;i++)
    {
          oldcs2[i]=getnumericalderiv(energy,pressure,i);
       //  myfile<<rhoB[i]/nsat<<" "<<energy[i]/(rhoB[i]*nsat)<<endl;
      //  myfile2<<rhoB[i]/nsat<<" "<<oldcs2[i]<<endl;
    }*/
    
     
    //get yQ information
    
   /*  double e0=31.7;
       double L=58.7;
       double K=-100;
       double J=32;
    ofstream myfile2("getYQ.txt");
    for(int i=0;i<limit;i++)
    {
        myfile2<<rhoB[i]/nsat<<" "<<getgausYQ(e0, L, K,rhoB[i])<<endl;
    }
*/
    
    /*Loop through different coefficients*/
    int count=0;
    double e0=27.5;
    double L=30;
    double K=-220;
    double J=-200;
    vector<vector<double>> all_cs2;
    vector<vector<double>> all_rhoB;
    ofstream myfile("all_cs2_eos3_35.txt");
   // ofstream myfile("all_energy_subtracted_sly2.txt");
    ofstream myfile2("coefficients_eos3_35.txt");
    int i1lim=8;
    int i2lim=10;
    int i3lim=8;
    int i4lim=9;
    bool causal=true;
    for(int i1=0;i1<i1lim;i1++)
    {
        for(int i2=0;i2<i2lim;i2++)
        {
           for(int i3=0;i3<i3lim;i3++)
           {
             //  std::cout<<" "<<rhoB[937]<<endl;
             //  cout<<K<<endl;
              for(int i4=0;i4<i4lim;i4++)
                {
                    double **cs2 = convert(rhoB,energy,e0,L,K,J,limit,0.5,0.15);
                    bool causal=isCausal(cs2[1],rhoB,0,4,5,limit);
                   // cout<<causal<<endl;
                  vector<double> single_cs2;
                    vector<double> single_rhoB;
                   if(causal)
                    {
                        for(int k=0;k<limit;k++)
                        {
                            single_cs2.push_back(cs2[1][k]);
                            single_rhoB.push_back(cs2[0][k]);
                        //    cout<<cs2[k]<<" ";
                        }
                        all_cs2.push_back(single_cs2);
                        all_rhoB.push_back(single_rhoB);
                        myfile2<<e0<<" "<<L<<" "<<K<<" "<<J<<endl;
                        cout<<e0<<" "<<L<<" "<<K<<" "<<J<<endl;
                       
                        count++;
                    }
                 //  else   cout<<e0<<" "<<L<<" "<<K<<" "<<J<<endl;
                    J=J+100;
                 
                }
                J=-250;
                K=K+40;
            //   cout<<e0<<" "<<L<<" "<<K<<" "<<J<<endl;
            }
            K=-220;
            L=L+10;
        }
        L=30;
        e0=e0+1;
        cout<<e0<<" "<<count<<endl;
    }
    //output cs2 information
   for(int i=0;i<count;i++)
    {
        for(int j=0;j<limit;j++)
        {
           // if(j=limit-1)
            //    myfile<<all_cs2[i][j]<<endl;
          //  else
               
         //   if(all_cs2[i][j]>2 || all_cs2[i][j]<-0.2) all_cs2[i][j]=0;
            myfile<<all_rhoB[i][j]<<" "<<all_cs2[i][j]<<endl;
          //  myfile<<all_cs2[i][j]<<endl;
          
        }
    }
    

//get single conversino information with specific symmetry energy coefficients
  /*    double e0=32.0;
      double L=47.46;
      double K=-115.13;
      double J=0;*/
     /*double e0=28;
      double L=30;
      double K=-100;
      double J=800;*/
    
 /*  double e0=30;
     double L=88;
     double K=0;
     double J=0;*/
  
    
/*  double **cs2=convert(rhoB,energy,e0,L,K,J,limit,0.5,0.15);
    ofstream myfile("yQ0_5_high.txt");
    int countt=0;
 //   ofstream myfile("nqmuq_eos3_01.txt");
    for(int i=0;i<limit;i++)
    {
      
            myfile<<cs2[0][i]/nsat<<" "<<cs2[1][i]<<endl;
    }*/

    
}

