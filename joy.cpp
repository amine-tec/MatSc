#include "../matplotlibcpp/matplotlibcpp.h"
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <vector>
#include <ctime>
#include <random>
#include <chrono>

using namespace std;
namespace plt = matplotlibcpp;

int main () {

   const double density =  7.13; //  g/cm**3
   const double atomicWeight = 65.38; //  u    
   const int Z = 30; 
         double Ei = 10.0; // 1.0 eV
         double J =  9.76*Z+ 58.52/pow(Z, 0.19); // mean ionization potential eV 
         double k =   0.734*pow(Z,0.037);      // parameter k                         
         std::vector<double> SP(1000); //eV/ang
         std::vector<double> E(SP.size());

         cout<<"Mean ionization potential(eV) = "<<J<<"and parameter K = "<<k<<"\n";

         for(size_t i=0; i<SP.size();i++){
         E[i]=Ei+i*10;
         SP[i] = 785.0*((density*Z)/(atomicWeight*E[i]))*log((1.166*(E[i]+k*J)/J));                        
         cout<<E[i]<<" eV "<<SP[i]<<" eV/ang"<<"\n";
         E[i]=(Ei+i*10)/1000;
         }
   

         plt::xlabel("Energy KeV");      
         plt::ylabel("Stopping power eV/ang");
         plt::title("Electron stopping power of Zn from Joy's ");
         plt::semilogx(E,SP);
         plt::show();
 
    
return 0;
}
