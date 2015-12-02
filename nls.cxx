#include <iostream>
#include <cmath>

using namespace std; 

// f(y)
void kutta(double *f, const double*  const y, const double  eta){
  f[0]=y[1];
  f[1]=(eta-y[0]*y[0])*y[0];
}
   
 

int main(){

const double dx=0.001;
const double X = 100;
const int N = X/dx;

const double eta=0.5;

double y[2];
y[0] = 1e-5;
y[1] = sqrt(eta) * y[0];

double k1[2];
double k2[2];
double k3[2];
double ytemp[2];


for(int i=0;i<N;i++){

  kutta(k1, y ,eta);  // k1 = f(y_n)
  
  ytemp[0] = y[0] + 0.5 * dx * k1[0];
  ytemp[1] = y[1] + 0.5 * dx * k1[1];
  
  kutta(k2, ytemp, eta); // k2 = f( y_n + 0.5 * dx * k1)
    
  ytemp[0] = y[0] + dx*(-k1[0]+2*k2[0]);
  ytemp[1] = y[1] + dx*(-k1[1]+2*k2[1]);
  
  kutta(k3,ytemp, eta); // k3 = f(y_n - dx * k1 + 2*dx*k2);
  
  y[0] = y[0] + dx/6 * (k1[0]+4*k2[0]+k3[0]);
  y[1] = y[1] + dx/6 * (k1[1]+4*k2[1]+k3[1]);

   cout << (i+1)*dx <<'\t'<< y[0]<< "\t" << y[1] <<endl;
  
}


 return 0;
}