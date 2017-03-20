#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstdio>

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
void analytical(double beta, double xMin, double xMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda);
double residual(double x, double lambda, double beta);
double qcrCalc(Truss &truss);
