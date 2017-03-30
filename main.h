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
double residual(double x, double lambda);
double KTExact(double x);
double qcrCalc(Truss &truss);
// tools.cpp
double power(double base, int exponent);
