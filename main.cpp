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
    int resolAna = 1000;
    double xMin = -0.5;
    double xMax = 2.5;
    double qcr = 2*12049281; // Factor 2 because 2*q w.r.t. notations in statement
    double qef = 1.2*qcr;

    // ----- Analytical solution -----
    //double qcr = qcrCalc(truss);
    std::vector<double> xAna;
    std::vector<double> lambdaAna;
    // Computes the analytical solution
    analytical(truss, xMin, xMax, resolAna, xAna, lambdaAna, qef);
    // Writes the results in a .txt file
    if(writeData(xAna, lambdaAna, "analytical")){return EXIT_FAILURE;}

    // ----- Incremental method -----
    double qefIncr = qef;
    int resolutionIncr = 11;
    std::vector<double> xIncr(resolutionIncr);
    std::vector<double> lambdaIncr(resolutionIncr);
    // Performs the incremental method
    incremental(truss, qefIncr, resolutionIncr, xIncr, lambdaIncr);
    // Writes the results in a .txt file
    if(writeData(xIncr, lambdaIncr, "incremental")){return EXIT_FAILURE;}

    // ----- Newton-Raphson method -----
    double qefNR = qef;
    int resolutionNR = 11;
    double epsilonNR = 0.0001;
    std::vector<double> xNR(resolutionNR);
    std::vector<double> lambdaNR(resolutionNR);
    // Performs the incremental method
    newtonRaphson(truss, qefNR, resolutionNR, epsilonNR, false, xNR, lambdaNR);
    // Writes the results in a .txt file
    if(writeData(xNR, lambdaNR, "NR")){return EXIT_FAILURE;}

    // ----- Spherical Arc-Length method -----
    double qefAL = qef;
    int maxIteration = 5;
    int idealIteration = 4;
    double epsilonAL = 0.001;
    double phi = 1;
    std::vector<double> xAL;
    std::vector<double> lambdaAL;
    arcLength(truss, qefAL, phi, maxIteration, idealIteration,
        epsilonAL, false, xAL, lambdaAL);

    //arcLength(truss, xMaxAL, dLambda0, phi, maxIteration, idealIteration,
    //    epsilonAL, false, xAL, lambdaAL);
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
