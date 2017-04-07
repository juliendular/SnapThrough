#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstdio>

#define BETA 2.598076211 // important parameter, equal to 9/(2*sqrt(3))


struct Bar{
    int nodes[2];
    double A;
    double E;
    double l0;
};

struct Truss{
    std::vector<double> nodes;
    std::vector<Bar> bars;
    std::vector<std::vector<int> > incidence;
    std::vector<bool> lockings;
    int nbDof;
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
// generalTruss.cpp
Truss simpleTruss();
Truss threeBarTruss();
Bar createBar(int nodeOne, int nodeTwo, double nodeOneX, double nodeOneY,
    double nodeTwoX, double nodeTwoY, double E, double A);
double PVW(Truss &truss, double u, double F);
void PVW(Truss &truss, std::vector<double> &u, std::vector<double> &F,
    std::vector<double> &OOB);


// simpleTruss.cpp
void analytical(Truss &truss, double xMin, double xMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda);
void incremental(Truss &truss, double lambdaMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda);
void newtonRaphson(Truss &truss, double lambdaMax, int resolution, double epsilon, bool modified,
    std::vector<double> &x, std::vector<double> &lambda);
void arcLength(Truss &truss, double xMax, double dLambda0, double phi,
    int maxIteration, int idealIteration,
    double epsilon, bool normal,
    std::vector<double> &x, std::vector<double> &lambda);

void predictorAL(Truss &truss, double x, double lambda, double phi, double Dl,
    double *dlambda, double *dx);


double residual(double x, double lambda);
double KTExact(double x);
double KTFD(Truss &truss, double x);
int invert(double *matrix, double *invMatrix, int n);
double qcrCalc(Truss &truss);
// tools.cpp
double power(double base, int exponent);
