#include "main.h"

int main(int argc, char *argv[]){

    // ----- Read input data -----
    if(argc == 3){
        //if(readData(argv[1], truss, param)){return EXIT_FAILURE;}
    }
    if(argc != 2){
        std::cout << "Wrong input arguments. Should be:" << std::endl;
        //std::cout << "\t- name of input data file" << std::endl;
        std::cout << "\t- name of ouput result file" << std::endl;
        return EXIT_FAILURE;
    }

    // ----- PARAMETERS -----
    // Truss
    Truss truss;
    truss.A = 0.1; // Cross-section area [m^2]
    truss.a = 1; // Horizontal projection of one bar [m]
    truss.b = 1; // Vertical projection of one bar [m]
    truss.E = 210000000000; // Young's modulus [N/m^2]
    truss.l0 = sqrt(truss.a*truss.a + truss.b*truss.b); // Length of one bar [m]
    // Numerical parameters
    int resolAna = 100;
    double xMin = -0.5;
    double xMax = 2.5;





    // ----- Analytical solution -----
    double qcr = qcrCalc(truss);
    std::vector <double> xAna;
    std::vector <double> lambda;

    double beta = 9.0/(2.0*sqrt(3));
    analytical(beta, xMin, xMax, resolAna, xAna, lambda);

    // Write the results in a .txt file
    if(writeData(xAna, lambda, "analytical")){return EXIT_FAILURE;}

    // ----- Incremental solution -----




    return 0;
}
