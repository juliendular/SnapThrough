#include "main.h"

int main(int argc, char *argv[]){

    // ----- Checks input argument -----
    if(argc != 1){
        std::cout << "Unknown input arguments. No arguments are needed" << std::endl;
        return EXIT_FAILURE;
    }

    // ----- PARAMETERS -----
    // Truss definition
    Truss truss;
    truss.A = 0.1; // Cross-section area [m^2]
    truss.a = 1; // Horizontal projection of one bar [m]
    truss.b = 1; // Vertical projection of one bar [m]
    truss.E = 210000000000; // Young's modulus [N/m^2]
    truss.l0 = sqrt(power(truss.a,2) + power(truss.b,2)); // Length of one bar [m]
    // Numerical parameters
    int resolAna = 100;
    double xMin = -0.5;
    double xMax = 2.5;

    // ----- Analytical solution -----
    double qcr = qcrCalc(truss);
    std::vector<double> xAna;
    std::vector<double> lambdaAna;
    // Computes the analytical solution
    analytical(xMin, xMax, resolAna, xAna, lambdaAna);
    // Writes the results in a .txt file
    if(writeData(xAna, lambdaAna, "analytical")){return EXIT_FAILURE;}

    // ----- Incremental solution -----
    double lambdaMaxIncr = 1.0;
    int resolutionIncr = 5;
    std::vector<double> xIncr(resolutionIncr);
    std::vector<double> lambdaIncr(resolutionIncr);
    // Performs the incremental method
    incremental(lambdaMaxIncr, resolutionIncr, xIncr, lambdaIncr);
    // Writes the results in a .txt file
    if(writeData(xIncr, lambdaIncr, "incremental")){return EXIT_FAILURE;}


    return 0;
}
