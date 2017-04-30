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
    int resolAna = 10000;
    double xMin = -0.5;
    double xMax = 2.5;
    double qcr = 2*12049281.52; // Factor 2 because 2*q w.r.t. notations in statement
    double qef = qcr;
/*
    // ----- Analytical solution -----
    std::vector<double> xAna;
    std::vector<double> lambdaAna;
    // Computes the analytical solution
    analytical(truss, xMin, xMax, resolAna, xAna, lambdaAna, qcr);
    // Writes the results in a .txt file
    if(writeData(xAna, lambdaAna, "analytical")){return EXIT_FAILURE;}

    // ----- Incremental method -----
    double qefIncr = qef;
    std::vector<int> resolutionIncr(3);
    resolutionIncr[0]=6; resolutionIncr[1]=11; resolutionIncr[2]=21;
    // Solving for some resolution values
    for(int i=0 ; i < resolutionIncr.size() ; i++){
        // Declaration
        std::vector<double> xIncr(resolutionIncr[i]);
        std::vector<double> lambdaIncr(resolutionIncr[i]);
        // Performs the incremental method
        incremental(truss, qefIncr, resolutionIncr[i], xIncr, lambdaIncr);
        // Determines the residuals
        std::vector<double> residualsIncr(resolutionIncr[i]);
        residuals(truss, xIncr, lambdaIncr, qef, residualsIncr);
        // Writes the results in a .txt file
        std::string file1 = "incremental" + std::to_string(i);
        std::string file2 = "incrementalResidual" + std::to_string(i);
        if(writeData(xIncr, lambdaIncr, file1)){return EXIT_FAILURE;}
        if(writeData(residualsIncr, file2)){return EXIT_FAILURE;}
    }
    // Above qcr
    int resolution = 21;
    std::vector<double> xIncrA(resolution);
    std::vector<double> lambdaIncrA(resolution);
    incremental(truss, 1.5*qefIncr, resolution, xIncrA, lambdaIncrA);
    if(writeData(xIncrA, lambdaIncrA, "incrementalAbove")){return EXIT_FAILURE;}

    // ----- Newton-Raphson method -----
    double qefNR = qef;
    std::vector<int> resolutionNR(3);
    resolutionNR[0]=6; resolutionNR[1]=11; resolutionNR[2]=21;
    double epsilonNR = 0.001;
    for(int i=0 ; i<resolutionNR.size() ; i++){
        std::vector<double> xNR(resolutionNR[i]);
        std::vector<double> lambdaNR(resolutionNR[i]);
        // Performs the incremental method
        newtonRaphson(truss, 0.0, qefNR, resolutionNR[i], epsilonNR, false, false, xNR, lambdaNR);
        // Determines the residuals
        std::vector<double> residualsNR(resolutionNR[i]);
        residuals(truss, xNR, lambdaNR, qefNR, residualsNR);
        // Writes the results in a .txt file
        std::string file1 = "NR" + std::to_string(i);
        std::string file2 = "NRResidual" + std::to_string(i);
        if(writeData(xNR, lambdaNR, file1)){return EXIT_FAILURE;}
        if(writeData(residualsNR, file2)){return EXIT_FAILURE;}
    }

    // ----- Displacement for a specific loading -----
    std::vector<double> x;
    std::vector<double> load;
    std::vector<double> stresses;
    int partialResolution = 31;
    std::vector<double> qefComposite(2); qefComposite[0] = qcr; qefComposite[1] = 2*qcr;
    for(int i=0 ; i<qefComposite.size() ; i++){
        // Definition on the piecewise linear loading
        std::vector<double> controlPoints(5);
        controlPoints[0] = 0.0; controlPoints[1] = qefComposite[i]; controlPoints[2] = 0.0;
        controlPoints[3] = -qefComposite[i]; controlPoints[4] = 0.0;
        // Numerical solver
        displacement(truss, partialResolution, epsilonNR, controlPoints, x, load, stresses);
        // Writes the results
        std::string file1 = "compositeLoading" + std::to_string(i);
        std::string file2 = "stresses" + std::to_string(i);
        if(writeData(x, load, file1)){return EXIT_FAILURE;}
        if(writeData(stresses, file2)){return EXIT_FAILURE;}
    }

    // ----- Spherical Arc-Length method -----
    double qefAL = 2.0*qcr;
    int maxIteration = 3;
    int idealIteration = 1;
    double epsilonAL = 0.001;
    double phi = 0.1/qefAL;
    double dLambdaInit = 0.1;
    std::vector<double> xAL;
    std::vector<double> lambdaAL;
    arcLength(truss, qefAL, phi, dLambdaInit, maxIteration, idealIteration,
        epsilonAL, 2, xAL, lambdaAL);
    // Writes the results in a .txt file
    if(writeData(xAL, lambdaAL, "AL")){return EXIT_FAILURE;}
*/
    // ----- Three bar truss -----
    Truss truss3 = threeBarTruss();
    std::vector<double> qef3(truss3.nbDof);
    qef3[0] = 0.0;
    qef3[1] = 2*qef;
    int maxIteration3 = 5;
    int idealIteration3 = 1;
    double epsilon3 = 0.01;
    double phi3 = 0.1/qef;
    double dLambdaInit3 = 0.02;
    std::vector<std::vector<double> > p3;
    std::vector<double> lambda3;
    // Arc-Length
    arcLength(truss3, qef3, phi3, dLambdaInit3,
        maxIteration3, idealIteration3, epsilon3, 0, p3, lambda3);
    // Write results
    if(writeData(p3, lambda3, "AL3")){return EXIT_FAILURE;}

/*
    // ----- Five/six bar truss -----
    Truss truss5 = fiveBarTruss();
    std::vector<double> qef5(truss5.nbDof);
    qef5[1] = 0.0;
    qef5[0] = 20*qef;
    int maxIteration5 = 5;
    int idealIteration5 = 1;
    double epsilon5 = 0.01;
    double phi5 = 0.1/qef;
    double dLambdaInit5 = 0.02;
    std::vector<std::vector<double> > p5;
    std::vector<double> lambda5;
    // Arc-Length
    arcLength(truss5, qef5, phi5, dLambdaInit5,
        maxIteration5, idealIteration5, epsilon5, 0, p5, lambda5);
    // Write results
    if(writeData(p5, lambda5, "AL5")){return EXIT_FAILURE;}

    Truss truss7 = sevenBarTruss();
    std::vector<double> qef7(truss7.nbDof);
    qef7[1] = 0.0;
    qef7[0] = 20*qef;
    int maxIteration7 = 5;
    int idealIteration7 = 1;
    double epsilon7 = 0.01;
    double phi7 = 0.1/qef;
    double dLambdaInit7 = 0.02;
    std::vector<std::vector<double> > p7;
    std::vector<double> lambda7;
    // Arc-Length
    arcLength(truss7, qef7, phi7, dLambdaInit7,
        maxIteration7, idealIteration7, epsilon7, 0, p7, lambda7);
    // Write results
    if(writeData(p7, lambda7, "AL7")){return EXIT_FAILURE;}

*/


    return 0;
}
