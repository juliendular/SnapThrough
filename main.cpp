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
/*
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
*//*
    // ----- Spherical Arc-Length method -----
    double qefAL = qef;
    int maxIteration = 5;
    int idealIteration = 4;
    double epsilonAL = 0.01;
    double phi = 0.3;
    double dLambdaInit = 0.07;
    std::vector<double> xAL;
    std::vector<double> lambdaAL;
    arcLength(truss, qefAL, phi, dLambdaInit, maxIteration, idealIteration,
        epsilonAL, false, xAL, lambdaAL);

    //arcLength(truss, 3, 0, phi, maxIteration, idealIteration,
    //    epsilonAL, false, xAL, lambdaAL);
    // Writes the results in a .txt file
    if(writeData(xAL, lambdaAL, "AL")){return EXIT_FAILURE;}
*/
    // Three bar truss

    //double qef = 2*12049281;

    Truss truss3 = threeBarTruss();

//*
    std::vector<double> qef3(truss3.nbDof);
    qef3[0] = 0;
    qef3[1] = 1.2*qef;
    int maxIteration3 = 3;
    int idealIteration3 = 2;
    double epsilon3 = 0.01;
    double phi3 = 0;
    double dLambdaInit3 = 0.04;
    std::vector<std::vector<double> > p3;
    std::vector<double> lambda3;
    // Arc-Length
    arcLength(truss3, qef3, phi3, dLambdaInit3,
        maxIteration3, idealIteration3, epsilon3, false, p3, lambda3);
    // Write results
    if(writeData(p3, lambda3, "AL3")){return EXIT_FAILURE;}
//*/
    // Test...
    std::vector<double> u(2);
    std::vector<double> F(2);
    std::vector<double> OOB(2,0.0);
    u[0] = 0.5; u[1] = 0.6;
    F[0] = 2e8; F[1] = 3e8;
    PVW(truss3, u, F, OOB);

    std::cout << "OOB force: ";
    for(int i=0 ; i<OOB.size() ; i++){
        std::cout << OOB[i] << " ";
    }
     std::cout << std::endl;

    return 0;
}
