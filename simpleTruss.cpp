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
void newtonRaphson(Truss &truss, double q0, double q1, int resolution, double epsilon,
    bool modified, bool start, std::vector<double> &x, std::vector<double> &lambda){
    // Charge increment
    double dlambda = 1/((double) resolution - 1);
    // First step: the equilibrium position for zero displacement/force
    if(!start)// If no specific starting point is given (otherwise x[0] should be set!)
        x[0] = 0.0;
    lambda[0] = 0.0; // Always starts from 0 and goes to 1
    // Loop until lambda = 1
    for(int it = 1 ; it < resolution ; it++){
        // Updates the force vector
        lambda[it] = lambda[it-1] + dlambda;
        // Computes the inverse tangent matrix
        double KT = KTFD(truss, x[it-1]);
        double KTinv;
        invert(&KT, &KTinv, 1);
        // Current load (linear interpolation)
        double load = (1-lambda[it])*q0 + lambda[it]*q1;
        // Predictor
        double xCurr = x[it-1] - KTinv * PVW(truss, x[it-1], load);
        // Corrector until epsilon precision is reached
        double currentResidual = PVW(truss, xCurr, load);
        int itCorr;
        for(itCorr=0 ; currentResidual > epsilon*abs(q0+q1) ||
              -currentResidual > epsilon*abs(q0+q1) ; itCorr++){
            // Computes the new inverse tangent matrix if not modified NR
            if(!modified){
                KT = KTFD(truss, xCurr);
                invert(&KT, &KTinv, 1);
            }
            // Corrects the displacement
            xCurr = xCurr - KTinv * currentResidual;
            currentResidual = PVW(truss, xCurr, load);
        }
        // Updates the displacement vector
        x[it] = xCurr;
    }
}

/* Displacement and stresses for piecewise linear loading for the 1 DOF case */
void displacement(Truss &truss, int partialResolution, double epsilon,
    std::vector<double> &controlPoints, std::vector<double> &x,
    std::vector<double> &load, std::vector<double> &stresses){
    // Starts from the equilibrium point
    x.resize(0); x.push_back(0.0);
    load.resize(0); load.push_back(0.0);
    stresses.resize(0); stresses.push_back(0.0);
    // Control points are assumed to be equidistant in time
    for(int i=0 ; i<controlPoints.size()-1 ; i++){
        // Declaration
        std::vector<double> partialX(partialResolution);
        std::vector<double> partialLambda(partialResolution);
        partialX[0] = x[x.size()-1];
        // Newton-Raphson method (with given starting point)
        newtonRaphson(truss, controlPoints[i], controlPoints[i+1], partialResolution, epsilon, false, true, partialX, partialLambda);
        // Insertion of the partial displacement
        x.insert(x.end(), partialX.begin()+1, partialX.end());
        // Determination of the partial load (not a direct output of NR method to remain consistent with other methods...)
        // and computation of the stresses
        for(int j=1 ; j<partialLambda.size() ; j++){
            load.push_back((1-partialLambda[j])*controlPoints[i] + partialLambda[j]*controlPoints[i+1]);
            stresses.push_back(stress(truss, partialX[j]));
        }
    }
}

/* Computation of the stresses for the simple 1 DOF truss */
double stress(Truss &truss, double x){
    //              E      / 2.0 * x * (x - 2*          b        )/         l0^2;
    return truss.bars[0].E / 2.0 * x * (x - 2*abs(truss.nodes[3]))/power(truss.bars[0].l0,2);
}



/* Calls the general arc-length method */
void arcLength(Truss &truss, double qef, double phi, double dLambdaInit, int maxIteration,
    int idealIteration, double epsilon, bool normal,
    std::vector<double> &x, std::vector<double> &lambda){
    // Declares one-element vectors
    std::vector<double> qefVec(1, qef);
    std::vector<std::vector<double> > p;
    // Calls the general function
    arcLength(truss, qefVec, phi, dLambdaInit, maxIteration, idealIteration, epsilon, normal,
        p, lambda);
    // Extract the information from p to x
    x.resize(p.size());
    for(int i=0 ; i<p.size() ; i++)
        x[i] = p[i][0];
    // Return
    return;
}

/* Implements the arc-length method */
void arcLength1(Truss &truss, double xMax, double dLambda0, double phi,
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
double residual(Truss &truss, double x, double lambda, double qef){
    double beta = truss.bars[0].E*truss.bars[0].A/(2.0*power(truss.bars[0].l0,3)*qef);
    return beta * x * (x-1) * (x-2) - lambda;
}

/* Computes the residuals for a list of load-displacement pairs */
void residuals(Truss &truss, std::vector<double> &x, std::vector<double> &lambda,
    double qef, std::vector<double> &res){
    // Resizes the residuals vector
    res.resize(x.size());
    // Computes the residual for each load/displacement pair
    for(int i=0 ; i<x.size() ; i++){
        res[i] = PVW(truss, x[i], lambda[i]*qef);
    }
}

/* Computes the analytical expression of the tangent stiffness matrix
IN: x: displacement ar which the tangent stiffness matrix has to be evaluated
OUT: the exact value of the tangent stifness matrix (scalar in this 1D case)
*/
double KTExact(Truss &truss, double x, double qef){
    double beta = truss.bars[0].E*truss.bars[0].A/(2.0*power(truss.bars[0].l0,3)*qef);
    return  beta * (3*x*x - 6*x + 2);
}

/* Computes the tangent stiffness matrix with a finite difference method */
double KTFD(Truss &truss, double x){
    double dx = 0.001;
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
