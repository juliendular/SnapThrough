#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstdio>

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
int writeData(std::vector<double> &x, std::string fileName);
int writeData(std::vector<double> &x, std::vector<double> &y, std::string filename);
int writeData(std::vector<std::vector<double> > &x, std::vector<double> &y, std::string fileName);
// generalTruss.cpp
Truss simpleTruss();
Truss threeBarTruss();
Bar createBar(int nodeOne, int nodeTwo, double nodeOneX, double nodeOneY,
    double nodeTwoX, double nodeTwoY, double E, double A);
double PVW(Truss &truss, double u, double F);
void PVW(Truss &truss, std::vector<double> &u, std::vector<double> &F, double lambda,
    std::vector<double> &OOB);
void PVW(Truss &truss, std::vector<double> &u, std::vector<double> &F,
    std::vector<double> &OOB);
// genaralSolver.cpp
void arcLength(Truss &truss, std::vector<double> qef, double phi, double dLambdaInit,
    int maxIteration, int idealIteration, double epsilon, bool normal,
    std::vector<std::vector<double> > &p, std::vector<double> &lambda);
void predictorAL(Truss &truss, std::vector<double> &p, double lambda, double phi,
    std::vector<double> &qef, double Dl, std::vector<double> &dp, double &dLambda);
void KTFD(Truss &truss, std::vector<double> &p, std::vector<std::vector<double> > &KT);
// simpleTruss.cpp
void analytical(Truss &truss, double xMin, double xMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda, double qef);
void incremental(Truss &truss, double qefIncr, int resolution,
    std::vector<double> &x, std::vector<double> &lambda);
void newtonRaphson(Truss &truss, double q0, double q1, int resolution, double epsilon, bool modified,
    bool start, std::vector<double> &x, std::vector<double> &lambda);
void displacement(Truss &truss, int partialResolution, double epsilon,
    std::vector<double> &controlPoints, std::vector<double> &x,
    std::vector<double> &load, std::vector<double> &stresses);
double stress(Truss &truss, double x);
void arcLength1(Truss &truss, double xMax, double dLambda0, double phi,
    int maxIteration, int idealIteration,
    double epsilon, bool normal,
    std::vector<double> &x, std::vector<double> &lambda);
void arcLength(Truss &truss, double qef, double phi, double dLambdaInit, int maxIteration,
    int idealIteration, double epsilon, bool normal,
    std::vector<double> &x, std::vector<double> &lambda);
void predictorAL(Truss &truss, double x, double lambda, double phi, double Dl,
    double *dlambda, double *dx);
double residual(Truss &truss, double x, double lambda, double qef);
void residuals(Truss &truss, std::vector<double> &x, std::vector<double> &lambda,
    double qef, std::vector<double> &res);
double KTExact(Truss &truss, double x, double qef);
double KTFD(Truss &truss, double x);
int invert(double *matrix, double *invMatrix, int n);
double qcrCalc(Truss &truss);
// tools.cpp
double power(double base, int exponent);
void vvPlus(std::vector<double> &v1, std::vector<double> &v2, std::vector<double> &sum);
void sv(double s, std::vector<double> &v, std::vector<double> &p);
double vv(std::vector<double> &v1, std::vector<double> &v2);
void mv(std::vector<std::vector<double> > &m, std::vector<double> &v,
    std::vector<double> &p);
int inv(std::vector<std::vector<double> > &m, std::vector<std::vector<double> > &invM);
bool posDefinite(std::vector<std::vector<double> > &m);
