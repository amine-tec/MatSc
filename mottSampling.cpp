 
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include <chrono>

using namespace std;

default_random_engine dre (chrono::steady_clock::now().time_since_epoch().count());     // provide seed
double random (double lim)
{
    uniform_real_distribution<double> uid {0,lim};   // help dre to generate nos from 0 to lim (lim included);
    return uid(dre);    // pass dre as an argument to uid to generate the random no
}


double getCosTheta(double angle[], double CuDeFun[]) // sampling the Mott Cross section
{

double min=0;
double thetaScatter=0;
double Random;

Random = random(1);

//std::srand(std::time(nullptr)); // use current time as seed for random generator

// Rand = std::rand();

// Rand=Rand/(RAND_MAX);


min=abs(Random-CuDeFun[0]);

int j=0;

for (int i=0; i<606; i++)

    {

  if (abs(Random-CuDeFun[i]) < min)
      
      {min = abs(Random-CuDeFun[i]);
          
        ++j;
       }

     }   

int lowerLimit;

if ((Random-CuDeFun[j])<=0)

{

//lowerLimit=j-1;

thetaScatter= ((angle[j]-angle[j-1])*Random-(angle[j])*CuDeFun[j-1]+(angle[j-1])*CuDeFun[j])/(CuDeFun[j]-CuDeFun[j-1]);

}

//lowerLimit=j;

else thetaScatter=thetaScatter= ((angle[j+1]-angle[j])*Random-(angle[j+1])*CuDeFun[j]+(angle[j])*CuDeFun[j+1])/(CuDeFun[j+1]-CuDeFun[j]);
 

// cout<<"Sampled Theta ="<<thetaScatter<<".\n";


return thetaScatter;
}

int main()
{

double ang[606];
double CDF[606];

ifstream infile;

infile.open("CDF.txt");

int num = 0;

while(!infile.eof()){
          
                     infile >> ang[num];
                     infile >> CDF[num];
                     ++num;

                     }

infile.close(); 


ofstream out_file("mott.txt");

  if (out_file.is_open()) {
for (int i=0; i<1000000 ; ++i)

                     {

cout<<""<<getCosTheta(ang,CDF)<<endl;

// plotAng[i] = getCosTheta(ang,CDF);

  out_file <<getCosTheta(ang,CDF)<<endl;
                     }

                           }

    out_file.close();

return 0;

}
