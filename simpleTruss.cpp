#include "main.h"

/* Fills the x and lambda vector with the analytical solution
IN: - xMin: starting displacement
    - xMax: final displacement
    - resolution: number of pairs (x, lambda) to be calculated. The x values
    are equidistant values between xMin and xMax
    - x: vector of displacement to be filled
    - lambda: vector of forces to be filled
*/
void analytical(Truss &truss, double xMin, double xMax, int resolution,
    std::vector<double> &x, std::vector<double> &lambda, double qef){

    double beta = truss.bars[0].E*truss.bars[0].A/(2.0*power(truss.bars[0].l0,3)*qef);
    double b = -truss.nodes[3];

    double dx = (xMax-xMin) / ((double)resolution - 1);
    double xCurrent = xMin;

    for(int i=0 ; i<resolution ; i++){
        x.push_back(xCurrent);
        // The factor 2* holds for the whole truss !!
        lambda.push_back(2 * beta * xCurrent * (xCurrent-b) * (xCurrent-2*b));
        xCurrent += dx;
    }

}

/* Implements the incremental method for a 1 DOF truss */
void incremental(Truss &truss, double qef, int resolution,
    std::vector<double> &x, std::vector<double> &lambda){
    // Charge increment
    double dlambda = 1/((double) resolution - 1);
    // First step: the equilibrium position for zero displacement/force
    x[0] = 0.0;
    lambda[0] = 0.0;
    // Loop until lambda = 1
    for(int it = 1 ; it < resolution ; it++){
        // Updates the force vector
        lambda[it] = lambda[it-1] + dlambda;
        // Computes the inverse tangent matrix
        double KT = KTFD(truss, x[it-1]);
        double KTinv;
        invert(&KT, &KTinv, 1);
        // Updates the displacement vector
        x[it] = x[it-1] - KTinv * PVW(truss, x[it-1], lambda[it]*qef);
    }
    return;
}

/* Implements the Newton-Raphson method for a 1 DOF truss */
void newtonRaphson(Truss &truss, double qef, int resolution, double epsilon,
    bool modified, std::vector<double> &x, std::vector<double> &lambda){
    // Charge increment
    double dlambda = 1/((double) resolution - 1);
    // First step: the equilibrium position for zero displacement/force
    x[0] = 0.0;
    lambda[0] = 0.0;
    // Loop until lambda = 1
    for(int it = 1 ; it < resolution ; it++){
        // Updates the force vector
        lambda[it] = lambda[it-1] + dlambda;
        // Computes the inverse tangent matrix
        double KT = KTFD(truss, x[it-1]);
        double KTinv;
        invert(&KT, &KTinv, 1);
        // Predictor
        double xCurr = x[it-1] - KTinv * PVW(truss, x[it-1], lambda[it]*qef);
        // Corrector until epsilon precision is reached
        double currentResidual = PVW(truss, xCurr, lambda[it]*qef);
        while(currentResidual > epsilon*lambda[it]*qef ||
              -currentResidual > epsilon*lambda[it]*qef){
            // Computes the new inverse tangent matrix if not modified NR
            if(!modified){
                KT = KTFD(truss, xCurr);
                invert(&KT, &KTinv, 1);
            }
            // Corrects the displacement
            xCurr = xCurr - KTinv * currentResidual;
            currentResidual = PVW(truss, xCurr, lambda[it]*qef);
        }
        // Updates the displacement vector
        x[it] = xCurr;
    }
}

/* Calls the general arc-length method */
void arcLength(Truss &truss, double qef, double phi, int maxIteration,
    int idealIteration, double epsilon, bool normal,
    std::vector<double> &x, std::vector<double> &lambda){
    // Declares one-element vectors
    std::vector<double> qefVec(1, qef);
    std::vector<std::vector<double> > p;
    // Calls the general function
    arcLength(truss, qefVec, phi, maxIteration, idealIteration, epsilon, normal,
        p, lambda);
    // Extract the information from p to x
    x.resize(p.size());
    for(int i=0 ; i<p.size() ; i++)
        x[i] = p[i][0];
    // Return
    return;
}

/* Implements the arc-length method */
void arcLength(Truss &truss, double xMax, double dLambda0, double phi,
    int maxIteration, int idealIteration,
    double epsilon, bool normal,
    std::vector<double> &x, std::vector<double> &lambda){
    // Initialization
    x.push_back(0.0);
    lambda.push_back(0.0);
    // Initial arc-length
    double Dl = 0.1; // USE NR INSTEAD
    // Loop until xMax
    for(int it = 1; x[it-1] < xMax ;){
        std::cout << "start iteration " << it << " ... ";
        // --- Predictor ---
        double dlambdaCurr;
        double dxCurr;
        predictorAL(truss, x[it-1], lambda[it-1], phi, Dl, &dlambdaCurr, &dxCurr);
        // Initialization of the corrector
        double xCurr = x[it-1] + dxCurr;
        double lambdaCurr = lambda[it-1] + dlambdaCurr;
        double currentResidual = PVW(truss, xCurr, lambdaCurr);
        // Loop on the corrector until convergence
        int corrIt;
        for(corrIt=1 ; corrIt < maxIteration && (currentResidual>epsilon || -currentResidual>epsilon) ; corrIt++){
            // --- Corrector ---
            // Computes the inverse of the tangent stiffness matrix
            double KT = KTFD(truss, xCurr);
            double KTinv;
            invert(&KT, &KTinv, 1);
            // Computes dxBar and dxt
            double dxBar = -KTinv * PVW(truss, xCurr, lambdaCurr);
            double dxt = KTinv;
            // Computes the a's coefficients
            double a1, a2, a3;
            a1 = power(dxt,2) + power(phi,2);
            a2 = 2*(dxt*(dxCurr+dxBar) + dlambdaCurr*power(phi,2));
            a3 = power((dxCurr+dxBar), 2) + power(dlambdaCurr*phi,2) - power(Dl,2);
            // Deduces the possible values for dlambda/dx
            double rho = power(a2,2) - 4*a1*a3;
            double dlambda1 = (-a2 + sqrt(rho))/(2*a1);
            double dlambda2 = (-a2 - sqrt(rho))/(2*a1);
            double dx1 = dxBar + dlambda1 * dxt;
            double dx2 = dxBar + dlambda2 * dxt;
            // Determines the good ones (no need to divide by Dl^2)
            double cos1 = (dxCurr*(dxCurr+dx1) + power(phi,2)*dlambdaCurr*(dlambdaCurr+dlambda1));
            double cos2 = (dxCurr*(dxCurr+dx2) + power(phi,2)*dlambdaCurr*(dlambdaCurr+dlambda2));
            if(cos1 > cos2){
                dlambdaCurr = dlambdaCurr + dlambda1;
                dxCurr = dxCurr + dx1;
            }
            else{
                dlambdaCurr = dlambdaCurr + dlambda2;
                dxCurr = dxCurr + dx2;
            }
            // Updates force and displacement
            xCurr = x[it-1] + dxCurr;
            lambdaCurr = lambda[it-1] + dlambdaCurr;
            currentResidual = PVW(truss, xCurr, lambdaCurr);
        }
        std::cout << "end with " << corrIt << " iterations" << std::endl;
        if(corrIt == maxIteration){
            Dl = Dl/4.0;
        }
        else{
            // Update arc-length
            Dl = Dl * sqrt( (double)idealIteration / (double)corrIt );
            // Save displacement and force
            x.push_back(xCurr);
            lambda.push_back(lambdaCurr);
            it++;
        }
    }

}

/* Provides predictor values for the arc-length method */
void predictorAL(Truss &truss, double x, double lambda, double phi, double Dl,
    double *dlambda, double *dx){
    // Computes the inverse of the tangent stiffness matrix
    double KT = KTFD(truss, x);
    double KTinv;
    invert(&KT, &KTinv, 1);
    // Computes dxBar and dxt
    double dxBar = -KTinv * PVW(truss, x, lambda);
    double dxt = KTinv;
    // Evaluate the positive definiteness of the tangent stiffness matrix
    bool positiveDefinite = (KT > 0); // TO IMPROVE
    // Computes the load increase predictor
    double term1 = - dxBar * dxt;
    double rad = power(dxBar*dxt,2)
                 + (power(Dl,2) - power(dxBar,2)) * (power(dxt,2)+power(phi,2));
    double deno = power(dxt,2) + power(phi,2);
    if(positiveDefinite)
        *dlambda = (term1 + sqrt(rad)) / deno;
    else
        *dlambda = (term1 - sqrt(rad)) / deno;
    // Deduce the displacement increase predictor
    *dx = dxBar + (*dlambda)*dxt;
    // Return
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

/* Computes the tangent stiffness matrix with a finite difference method */
double KTFD(Truss &truss, double x){
    double dx = 0.001; // GOOD VALUE ???
    return (PVW(truss, x+dx, 0) - PVW(truss, x-dx, 0)) / (2*dx);
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


/* TO CHANGE !!!
// Computes the critical charge (in [N] !) of the truss given as a parameter
double qcrCalc(Truss &truss){
    return sqrt(3)*truss.A*truss.E*power(truss.b, 3) / (9.0*power(truss.l0,3));
}
*/
