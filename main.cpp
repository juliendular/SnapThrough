#include "main.h"

int main(int argc, char *argv[]){

    // ----- Read input data -----
    Truss truss;
    Parameters param;
    if(argc == 3){
        if(readData(argv[1], truss, param)){return EXIT_FAILURE;}
    }
    else{
        std::cout << "Wrong input arguments. Should be:" << std::endl;
        std::cout << "\t- name of input data file" << std::endl;
        std::cout << "\t- name of ouput result file" << std::endl;
        return EXIT_FAILURE;
    }

    // ----- Analytical solution -----
    double qcr = qcrCalc(truss);
    std::vector <double> xAna;
    std::vector <double> lambda;

    double beta = 9.0/(2.0*sqrt(3));
    analytical(beta, param.xMin, param.xMax, param.resolAna, xAna, lambda);

    // Write the data in a .txt file
    if(writeData(xAna, lambda, "analytical")){return EXIT_FAILURE;}

    // ----- Incremental solution -----




    return 0;
}
