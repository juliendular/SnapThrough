#include "main.h"

int main(int argc, char *argv[]){

    // Read input data from file
    Truss truss;
    Parameters parameters;
    if(argc == 3){
        if(readData(argv[1], truss, parameters)){
            std::cout << "Invalid input file" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else{
        std::cout << "Wrong input arguments. Should be:" << std::endl;
        std::cout << "\t- name of input data file" << std::endl;
        std::cout << "\t- name of ouput result file" << std::endl;
        return EXIT_FAILURE;
    }

    // Analytical solution
    std::cout << "Analytical solution" << std::endl;

    std::cout << truss.A << " " << truss.l0 << std::endl;
    std::cout << parameters.resolution << " " << parameters.tFinal << std::endl;





    // Incremental solution
    std::cout << "Incremental solution" << std::endl;


    return 0;
}
