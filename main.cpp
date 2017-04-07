#include "main.h"

int main(int argc, char *argv[]){

    // ----- Checks input argument -----
    if(argc != 1){
        std::cout << "Unknown input arguments. No arguments are needed" << std::endl;
        return EXIT_FAILURE;
    }

    // ----- Simple truss study -----
    Truss truss = simpleTruss();

    // ----- Parameters -----
    // Numerical parameters
    int resolAna = 100;
    double xMin = -0.5;
    double xMax = 2.5;

    // ----- Analytical solution -----
    //double qcr = qcrCalc(truss);
    std::vector<double> xAna;
    std::vector<double> lambdaAna;
    // Computes the analytical solution
    analytical(truss, xMin, xMax, resolAna, xAna, lambdaAna);
    // Writes the results in a .txt file
    if(writeData(xAna, lambdaAna, "analytical")){return EXIT_FAILURE;}

    // ----- Incremental method -----
    double lambdaMaxIncr = 1.5;
    int resolutionIncr = 10;
    std::vector<double> xIncr(resolutionIncr);
    std::vector<double> lambdaIncr(resolutionIncr);
    // Performs the incremental method
    incremental(truss, lambdaMaxIncr, resolutionIncr, xIncr, lambdaIncr);
    // Writes the results in a .txt file
    if(writeData(xIncr, lambdaIncr, "incremental")){return EXIT_FAILURE;}

    // ----- Newton-Raphson method -----

    double lambdaMaxNR = 1.5;
    int resolutionNR = 10;
    double epsilonNR = 0.01;
    std::vector<double> xNR(resolutionNR);
    std::vector<double> lambdaNR(resolutionNR);
    // Performs the incremental method
    newtonRaphson(truss, lambdaMaxNR, resolutionNR, epsilonNR, false, xNR, lambdaNR);
    // Writes the results in a .txt file
    if(writeData(xNR, lambdaNR, "NR")){return EXIT_FAILURE;}

    // ----- Spherical Arc-Length method -----
    double xMaxAL = 2.5;
    int maxIteration = 5;
    int idealIteration = 4;
    double dLambda0 = 0.1;
    double epsilonAL = 0.01;
    double phi = 1;
    std::vector<double> xAL;
    std::vector<double> lambdaAL;
    arcLength(truss, xMaxAL, dLambda0, phi, maxIteration, idealIteration,
        epsilonAL, false, xAL, lambdaAL);
    // Writes the results in a .txt file
    if(writeData(xAL, lambdaAL, "AL")){return EXIT_FAILURE;}





    // Test...
    std::vector<double> u(2);
    std::vector<double> F(2);
    std::vector<double> OOB(2,0.0);
    u[0] = 0; u[1] = 0.1;
    F[0] = 0; F[1] = 0;
    PVW(truss, u, F, OOB);

    std::cout << "OOB force: " << PVW(truss, 0.1, 0) << std::endl;

    return 0;
}
