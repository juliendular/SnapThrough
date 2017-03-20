#include "main.h"

void analytical(double beta, double xMin, double xMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda){

    double dx = (xMax-xMin) / (double)resolution;
    double xCurrent = xMin;

    for(int i=0 ; i<resolution ; i++){
        x.push_back(xCurrent);
        lambda.push_back(beta * xCurrent * (xCurrent-1) * (xCurrent-2));
        xCurrent += dx;
    }

}

double residual(double x, double lambda, double beta){
    return beta * x * (x-1) * (x-2) - lambda;
}

double qcrCalc(Truss &truss){
    return sqrt(3)*truss.A*truss.E*truss.b*truss.b*truss.b/(9.0*truss.l0*truss.l0*truss.l0);
}
/*
double beta(Truss &truss){

}
*/
