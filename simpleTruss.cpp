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
    // Charge increment
    double dlambda = lambdaMax/((double) resolution);
    // First step: the equilibrium position for zero displacement/force
    x[0] = 0.0;
    lambda[0] = 0.0;
    // Loop until lambdaMax
    for(int it = 1 ; it < resolution ; it++){
        // Updates the force vector
        lambda[it] = lambda[it-1] + dlambda;
        // Computes the inverse tangent matrix
        double KT = KTExact(x[it-1]);
        double KTinv;
        invert(&KT, &KTinv, 1);
        // Updates the displacement vector
        x[it] = x[it-1] - KTinv * residual(x[it-1], lambda[it]);
    }
    return;
}

/* Implements the Newton-Raphson method */
void newtonRaphson(double lambdaMax, int resolution, double epsilon, bool modified,
    std::vector<double> &x, std::vector<double> &lambda){
    // Charge increment
    double dlambda = lambdaMax/((double) resolution);
    // First step: the equilibrium position for zero displacement/force
    x[0] = 0.0;
    lambda[0] = 0.0;
    // Loop until lambdaMax
    for(int it = 1 ; it < resolution ; it++){
        // Updates the force vector
        lambda[it] = lambda[it-1] + dlambda;
        // Computes the inverse tangent matrix
        double KT = KTExact(x[it-1]);
        double KTinv;
        invert(&KT, &KTinv, 1);
        // Predictor
        double xCurr = x[it-1] - KTinv * residual(x[it-1], lambda[it]);
        // Corrector until epsilon precision is reached
        double currentResidual = residual(xCurr, lambda[it]);
        while(currentResidual > epsilon || -currentResidual > epsilon){
            // Computes the new inverse tangent matrix if not modified NR
            if(!modified){
                KT = KTExact(xCurr);
                invert(&KT, &KTinv, 1);
            }
            // Corrects the displacement
            xCurr = xCurr - KTinv * currentResidual;
            currentResidual = residual(xCurr, lambda[it]);
        }
        // Updates the displacemet vector
        x[it] = xCurr;
    }
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

/* Inverts a n*n matrix */
int invert(double *matrix, double *invMatrix, int n){
    if(n==1){
        if(*matrix == 0.0){
            std::cout << "Division by zero" << std::endl;
            return 1;
        }
        else{
            *invMatrix = 1.0/(*matrix);
            return 0;
        }
    }
    else if(n==2){
        std::cout << "2*2 matrix inversion not implemented" << std::endl;
        return 0; // TO IMPLEMENT
    }
    else{
        std::cout << "n*n (n>2) matrix inversion not implemented" << std::endl;
        return 0; // TO IMPLEMENT
    }
}



// Computes the critical charge (in [N] !) of the truss given as a parameter
double qcrCalc(Truss &truss){
    return sqrt(3)*truss.A*truss.E*power(truss.b, 3) / (9.0*power(truss.l0,3));
}
