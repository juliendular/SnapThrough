#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <ctime>

struct Truss{
    double A, a, b, l0, E;
};

struct Parameters{
    int resolution;
    double tFinal;
};

int readData(std::string fileName, Truss &truss, Parameters &parameters);
