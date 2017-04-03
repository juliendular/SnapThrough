#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstdio>

#define BETA 2.598076211 // important parameter, equal to 9/(2*sqrt(3))

struct Truss{
    double A, a, b, l0, E;
};

struct Parameters{
    int resolAna;
    double xMin;
    double xMax;
};

// readData.cpp
int readData(std::string fileName, Truss &truss, Parameters &parameters);
// writeData.cpp
int writeData(std::vector<double> &x, std::vector<double> &y, std::string filename);
// simpleTruss.cpp
void analytical(double xMin, double xMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda);
void incremental(double lambdaMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda);
void newtonRaphson(double lambdaMax, int resolution, double epsilon, bool modified,
    std::vector<double> &x, std::vector<double> &lambda);
void arcLength(double xMax, double dLambda0, double phi,
    int maxIteration, int idealIteration,
    double epsilon, bool normal,
    std::vector<double> &x, std::vector<double> &lambda);

void predictorAL(double x, double lambda, double phi, double Dl,
    double *dlambda, double *dx);


double residual(double x, double lambda);
double KTExact(double x);
int invert(double *matrix, double *invMatrix, int n);
double qcrCalc(Truss &truss);
// tools.cpp
double power(double base, int exponent);
