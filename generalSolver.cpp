#include "main.h"

/* Implements the Newton-Raphson method for a multidimensinal truss */
double initAL(Truss &truss, std::vector<double> &qef, double dLambda, double epsilon, double phi){
    // First step: the equilibrium position for zero displacement/force
    std::vector<double> p(truss.nbDof,0.0);
    double lambda = 0.0 + dLambda;
    // Computes the inverse of the tangent stiffness matrix
    std::vector<std::vector<double> > KT(truss.nbDof);
    std::vector<std::vector<double> > KTinv(truss.nbDof);
    for(int i=0 ; i<truss.nbDof ; i++){
        KT[i].resize(truss.nbDof);
        KTinv[i].resize(truss.nbDof);
    }
    KTFD(truss, p, KT);
    inv(KT, KTinv);
    // Predictor
    std::vector<double> OOB(truss.nbDof);
    std::vector<double> dp(truss.nbDof);
    PVW(truss, p, qef, lambda, OOB);
    mv(KTinv, OOB, dp);
    sv(-1, dp, dp);
    vvPlus(p, dp, p);
    PVW(truss, p, qef, lambda, OOB);
    double OOBeq = sqrt(vv(OOB, OOB)/vv(qef,qef))/lambda;
    // Corrector until epsilon precision is reached
    while(OOBeq > epsilon || -OOBeq > epsilon){
        // Does not compute the new inverse tangent matrix (modified NR)
        // Corrects the displacement
        mv(KTinv, OOB, dp);
        sv(-1, dp, dp);
        vvPlus(p, dp, p);
        PVW(truss, p, qef, lambda, OOB);
        OOBeq = sqrt(vv(OOB, OOB)/vv(qef,qef));
    }
    // Deduce the arc-length
    return sqrt(vv(p,p) + power(phi*dLambda,2)*vv(qef,qef));
}

/* Implements the multidimensional arc-length method */
void arcLength(Truss &truss, std::vector<double> qef, double phi, double dLambdaInit,
    int maxIteration, int idealIteration, double epsilon, int normal,
    std::vector<std::vector<double> > &p, std::vector<double> &lambda){
    // Initialization
    int nbDof = truss.nbDof;
    std::vector<double> pInit(nbDof, 0.0);
    p.push_back(pInit);
    lambda.push_back(0.0);
    // Initial arc-length
    double Dl = initAL(truss, qef, dLambdaInit, epsilon, phi);
    // Loop until lambda = 1
    int it = 1;
    while(lambda[it-1]<1){
        // --- Predictor ---
        double dLambdap; // For the normal arc-length method
        double dLambda0;
        std::vector<double> dpp(nbDof); // For the normal arc-length method
        std::vector<double> dp0(nbDof);
        predictorAL(truss, p[it-1], lambda[it-1], phi, qef, Dl, dpp, dLambdap);
        dLambda0 = dLambdap;
        dp0 = dpp;
        // Initialization of the corrector
        std::vector<double> p0;
        vvPlus(p[it-1], dp0, p0);
        double lambda0 = lambda[it-1] + dLambda0;
        std::vector<double> OOB(nbDof);
        PVW(truss, p0, qef, lambda0, OOB);

        // FOR INFO
        std::vector<double> testsP(0); testsP.push_back(p0[0]);
        std::vector<double> testsL(0); testsL.push_back(lambda0);


        // Estimate the adimensional OOB force
        double OOBeq = sqrt(vv(OOB, OOB)/vv(qef,qef));
        // Loop on the corrector until convergence
        int corrIt;
        for(corrIt=1 ; corrIt < maxIteration && (OOBeq>epsilon || -OOBeq>epsilon) ; corrIt++){
            // --- Corrector ---
            // Computes the inverse of the tangent stiffness matrix
            std::vector<std::vector<double> > KT(truss.nbDof);
            std::vector<std::vector<double> > KTinv(truss.nbDof);
            for(int i=0 ; i<truss.nbDof ; i++){
                KT[i].resize(truss.nbDof);
                KTinv[i].resize(truss.nbDof);
            }
            KTFD(truss, p0, KT);
            inv(KT, KTinv);
            // Computes dxBar and dxt
            std::vector<double> dpt(truss.nbDof);
            std::vector<double> dpBar(truss.nbDof);
            mv(KTinv, qef, dpt);
            mv(KTinv, OOB, dpBar);
            sv(-1.0, dpBar, dpBar);
            // Detemines the method
            if(normal==0){ // Spherical arc-length
                // Computes the a's coefficients
                double a1, a2, a3;
                std::vector<double> tmp(nbDof);
                vvPlus(dp0, dpBar, tmp);
                a1 = vv(dpt, dpt) + power(phi,2) * vv(qef, qef);
                a2 = 2*vv(dpt, tmp) + 2*dLambda0*power(phi,2)*vv(qef, qef);
                a3 = vv(tmp, tmp) - power(Dl,2) + power(dLambda0*phi, 2)*vv(qef, qef);
                // Deduces the possible values for dlambda/dx
                double rho = power(a2,2) - 4*a1*a3;
                double dLambda1 = (-a2 + sqrt(rho))/(2*a1);
                double dLambda2 = (-a2 - sqrt(rho))/(2*a1);
                std::vector<double> dp1(nbDof); // = dpBar + dlambda1 * dpt;
                std::vector<double> dp2(nbDof); // = dpBar + dlambda2 * dpt;
                sv(dLambda1, dpt, tmp);
                vvPlus(dpBar, tmp, dp1);
                sv(dLambda2, dpt, tmp);
                vvPlus(dpBar, tmp, dp2);
                // Determines the good ones (no need to divide by Dl^2)
                std::vector<double> dpNext1(nbDof);
                std::vector<double> dpNext2(nbDof);
                vvPlus(dp0, dp1, dpNext1);
                vvPlus(dp0, dp2, dpNext2);
                double cos1 = vv(dp0,dpNext1) + power(phi,2)*vv(qef,qef)*dLambda0*(dLambda0+dLambda1);
                double cos2 = vv(dp0,dpNext2) + power(phi,2)*vv(qef,qef)*dLambda0*(dLambda0+dLambda2);
                if(cos1 > cos2){
                    dLambda0 += dLambda1;
                    dp0 = dpNext1;
                }
                else{
                    dLambda0 += dLambda2;
                    dp0 = dpNext2;
                }
            }
            else if(normal == 1){ // Updated normal arc-length
                double dLambda1 = - vv(dp0, dpBar) / (vv(dp0,dpt) + dLambda0*power(phi,2)*vv(qef,qef));
                std::vector<double> dp1(nbDof); // = dpBar + dlambda1 * dpt;
                sv(dLambda1, dpt, dp1);
                vvPlus(dpBar, dp1, dp1);
                // Saving
                dLambda0 += dLambda1;
                vvPlus(dp0, dp1, dp0);
            }
            else{ // Normal arc-length
                double dLambda1 = - vv(dpp, dpBar) / (vv(dpp,dpt) + dLambdap*power(phi,2)*vv(qef,qef));
                std::vector<double> dp1(nbDof); // = dpBar + dlambda1 * dpt;
                sv(dLambda1, dpt, dp1);
                vvPlus(dpBar, dp1, dp1);
                // Saving
                dLambda0 += dLambda1;
                vvPlus(dp0, dp1, dp0);
            }
            // Updates force and displacement
            vvPlus(p[it-1], dp0, p0);
            lambda0 = lambda[it-1] + dLambda0;
            PVW(truss, p0, qef, lambda0, OOB);
            OOBeq = sqrt(vv(OOB, OOB)/vv(qef,qef));

            // FOR INFO
            testsP.push_back(p0[0]);
            testsL.push_back(lambda0);

        }

        // FOR INFO
        std::string name = "info" + std::to_string(it);// + "_" + std::to_string(corrIt);
        writeData(testsP, testsL, name);

        if(corrIt == maxIteration){
            Dl = Dl/2.0;
            std::cout << "Restart step " << it << std::endl;
        }
        else{
            std::cout << "End step " << it << " with " << corrIt << " corrector iterations." << std::endl;
            // Update arc-length
            Dl = Dl * sqrt( (double)idealIteration / (double)corrIt );
            // Save displacement and force
            p.push_back(p0);
            lambda.push_back(lambda0);
            it++;
        }
    }
}

/* Computes a predictor for the multidimensional arc-length method */
void predictorAL(Truss &truss, std::vector<double> &p, double lambda, double phi,
    std::vector<double> &qef, double Dl, std::vector<double> &dp, double &dLambda){
    // Computes the inverse of the tangent stiffness matrix
    std::vector<std::vector<double> > KT(truss.nbDof);
    std::vector<std::vector<double> > KTinv(truss.nbDof);
    for(int i=0 ; i<truss.nbDof ; i++){
        KT[i].resize(truss.nbDof);
        KTinv[i].resize(truss.nbDof);
    }
    KTFD(truss, p, KT);
    inv(KT, KTinv);

    // Deduces dpt (assumption residual is small)
    std::vector<double> dpt(truss.nbDof);
    std::vector<double> dpBar(truss.nbDof);
    std::vector<double> OOB(truss.nbDof);
    PVW(truss, p, qef, lambda, OOB);
    mv(KTinv, qef, dpt);
    mv(KTinv, OOB, dpBar);
    sv(-1.0, dpBar, dpBar);
    // Evaluate the positive definiteness of KT
    bool posDef = posDefinite(KT);
    // Computes the load increase
    double term1 = - vv(dpBar, dpt);
    double rad = vv(dpBar,dpBar)*vv(dpt,dpt)
                 + (power(Dl,2) - vv(dpBar,dpBar)) * (vv(dpt,dpt)+power(phi,2)*vv(qef,qef));
    double deno = vv(dpt,dpt) + power(phi,2)*vv(qef,qef);
    if(posDef){dLambda = (term1 + sqrt(rad)) / deno;}
    else{dLambda = (term1 - sqrt(rad)) / deno;}
    /*
    double rad = vv(dpt, dpt) + power(phi,2) * vv(qef, qef);
    if(posDef){dLambda = Dl/sqrt(rad);}
    else{dLambda = -Dl/sqrt(rad);}//*/
    // Deduce the displacement increase
    sv(dLambda, dpt, dp);
    vvPlus(dp, dpBar, dp);
    // Return
    return;
}

/* Computes the tangent stifness matrix */
void KTFD(Truss &truss, std::vector<double> &p, std::vector<std::vector<double> > &KT){
    // Step for central finite differenc
    double dp = 0.001;
    // Vector of zero forces
    std::vector<double> F(truss.nbDof, 0.0);
    // Stiffness for each direction
    for(int i=0 ; i<KT.size() ; i++){
        // Displacement
        std::vector<double> pP(truss.nbDof);
        std::vector<double> pM(truss.nbDof);
        pP = p;
        pM = p;
        pP[i] += dp/2.0;
        pM[i] += -dp/2.0;
        // OOB forces
        std::vector<double> OOB_P(truss.nbDof);
        std::vector<double> OOB_M(truss.nbDof);
        std::vector<double> F(truss.nbDof, 0.0);
        PVW(truss, pP, F, OOB_P);
        PVW(truss, pM, F, OOB_M);
        // Finite difference partial derivative
        std::vector<double> stiffness(truss.nbDof);
        sv(-1, OOB_M, OOB_M);
        vvPlus(OOB_P, OOB_M, stiffness);
        sv(1.0/dp, stiffness, stiffness);
        // Saving in KT
        KT[i] = stiffness;
    }
}
