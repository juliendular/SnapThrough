#include "main.h"

/* Fills the x and lambda vector with the analytical solution
IN: - xMin: starting displacement
    - xMax: final displacement
    - resolution: number of pairs (x, lambda) to be calculated. The x values
    are equidistant values between xMin and xMax
    - x: vector of displacement to be filled
    - lambda: vector of forces to be filled
*/
void analytical(double xMin, double xMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda){

    double dx = (xMax-xMin) / (double)resolution;
    double xCurrent = xMin;

    for(int i=0 ; i<resolution ; i++){
        x.push_back(xCurrent);
        lambda.push_back(BETA * xCurrent * (xCurrent-1) * (xCurrent-2));
        xCurrent += dx;
    }

}

/* Implements the incremental method */
void incremental(double lambdaMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda){
    // Charge incremental
    double dlambda = lambdaMax/((double) resolution);
    // First step: the equilibrium position for zero displacement/force
    x[0] = 0.0;
    lambda[0] = 0.0;
    // Loop until lambdaMax
    for(int it = 1 ; it < resolution ; it++){
        // Updates the force vector
        lambda[it] = lambda[it-1] + dlambda;
        // Updates the displacement vector
        x[it] = x[it-1] - 1/KTExact(x[it-1]) * residual(x[it-1], lambda[it]);
    }
    return;
}

/* Computes the residual force corresponding to the state (x, lambda)
IN: - x: displacement
    - lambda: force
OUT: the residual/out-of-balance force corresponding to the pair (x,lambda)
*/
double residual(double x, double lambda){
    return BETA * x * (x-1) * (x-2) - lambda;
}

/* Computes the analytical expression of the tangent stiffness matrix
IN: x: displacement ar which the tangent stiffness matrix has to be evaluated
OUT: the exact value of the tangent stifness matrix (scalar in this 1D case)
*/
double KTExact(double x){
    return BETA * (3*x*x - 6*x + 2);
}

// Computes the critical charge (in [N] !) of the truss given as a parameter
double qcrCalc(Truss &truss){
    return sqrt(3)*truss.A*truss.E*power(truss.b, 3) / (9.0*power(truss.l0,3));
}
