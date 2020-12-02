#include "matplotlibcpp.h"
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <vector>
#include <iterator>
#include <ctime>
#include <random>
#include <chrono>


using namespace std;
namespace plt = matplotlibcpp;

double ang[606];
double CDF[606];

default_random_engine dre (chrono::steady_clock::now().time_since_epoch().count());     
double random (double lim)
{
	uniform_real_distribution<double> uid {0,lim};   
	return uid(dre);    
}

double getCosTheta(const double angle[], const double CuDeFun[]) // sampling the Mott Cross section
{
	double min=0;
	double thetaScatter=0;
	double Rand=0;

	Rand=random(1);
	min=(Rand-CuDeFun[0]);

	int j=0;
    for(int i=0; i<606; i++){
		if (((Rand-CuDeFun[i]) < min)&& ((Rand-CuDeFun[i])>=0)){
		min = abs(Rand-CuDeFun[i]);     
		++j;
		}
	}   
	
	thetaScatter=thetaScatter= ((angle[j+1]-angle[j])*Rand-(angle[j+1])*CuDeFun[j]+(angle[j])*CuDeFun[j+1])/(CuDeFun[j+1]-CuDeFun[j]);

return cos(thetaScatter*M_PI/180);

}

void spin(double &mu_x, double &mu_y, double &mu_z)
{
    double costheta = getCosTheta(ang,CDF);
    double phi = 2 * M_PI * random(1);
    double sintheta = sqrt(1.0 - costheta * costheta); // sin(theta)
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    if (mu_z == 1.0) {
        mu_x = sintheta * cosphi;
        mu_y = sintheta * sinphi;
        mu_z = costheta;
    }
    else if (mu_z == -1.0) {
        mu_x = sintheta * cosphi;
        mu_y = -sintheta * sinphi;
        mu_z = -costheta;
    }
    else {
        double denom = sqrt(1.0 - mu_z * mu_z);
        double muzcosphi = mu_z * cosphi;
        double ux = sintheta * (mu_x * muzcosphi - mu_y * sinphi) / denom + mu_x * costheta;
        double uy = sintheta * (mu_y * muzcosphi + mu_x * sinphi) / denom + mu_y * costheta;
        double uz = -denom * sintheta * cosphi + mu_z * costheta;
        mu_x = ux, mu_y = uy, mu_z = uz;
    }

}

int main(){ 

    ifstream infile;
    infile.open("CDF.txt");    
    int num = 0;
    double mu_x,mu_y=0;
    double mu_z=1;
    while (!infile.eof()){      
          infile >> ang[num];
          infile >> CDF[num];
          ++num;
          }
    infile.close();

    // compute diffuse reflectance and transmittance
    const int nphotons = 100000; 
    double scale = 1.0 / nphotons; 
    double sigma_a = 0, sigma_s = 5.237*pow(5.29e-9,2), sigma_t = sigma_a + sigma_s; 
    double d = 1000;
    const int m = 10;
    double Rd = 0, Tt = 0;
    for (int n = 0; n < nphotons; ++n) { 
        double w = 1; 
        double x = 0, y = 0, z = 0, mux = 0, muy = 0, muz = 1, E = 1; 
//		cout<<0.0<<", "<<0.0<<", "<<0.0<<""<<"\n"; 
        std::vector<double> energyBack; 
		double Na = 6.023E+23, Z = 30, density = 7.13, atomicWeight = 65.38 , k = 0.734*pow(Z,0.037), J = (9.76*Z+ 58.52/pow(Z, 0.19))*1.0E-03;
        while (1) { 
            double s = -log(random(1)) /((Na*density/atomicWeight)*sigma_t); 
            double distToBoundary = 0; 
            if (muz > 0) distToBoundary = (d - z) / muz; 
            else if (muz < 0) distToBoundary = -z / muz; 
            if (s > distToBoundary) { // if(E<50) break; if energy bellow 50 eV don't enter loop
                if (muz > 0) Tt += w; else Rd += w;           cout<<x*1.0E+08<<", "<<y*1.0E+08<<", "<<E<<"\n";//  cout<<E<<"\n";//energyBack.push_back(E); 
                break; 
            } 
            // move
            x += s * mux; 
            y += s * muy; 
            z += s * muz; 
            E =  E-s*7.85E+04*(density*Z)/(atomicWeight*E)*log((1.166*(E+k*J))/J);  if(E<50.0E-03) break;
           
            // absorb if energy < 50.0 eV
/*
            double dw = sigma_a / sigma_t; 
            w -= dw; w = std::min(0.0, w); 
            if (w < 0.001) { // russian roulette test 
                if (random(1) > 1.0 / m) break; 
                else w *= m; 
            }
*/ 
            // scatter
            spin(mux, muy, muz); 
          //  cout<<random(1)<<".\n";
            
 //           cout<<x*1.0E+08<<", "<<y*1.0E+08<<", "<<z*1.0E+08<<"\n";  
//         X.push_back(x); Y.push_back(y); Zi.push_back(z); 
//            plt::fig(;
//            plt::axes(projection="3d");
            

                            
        }
/*    
     plt::figure(111);
     plt::plot_surface(X,Y,Z);
     plt::show();
     X.clear();
     Y.clear();
     Zi.clear();
*/
    } 

// cout<<"Transmittance = "<<Tt*scale<<" Diffuse reflectance = "<<Rd*scale<<" eV"<<"\n";  



return 0;

}
