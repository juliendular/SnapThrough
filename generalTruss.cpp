#include "main.h"

Truss simpleTruss(){
    Truss truss;
    // Field declaration
    std::vector<double> nodes;
    std::vector<Bar> bars;
    std::vector<std::vector<int> > incidence;
    std::vector<bool> lockings;
    // Field declaration
    int nbNodes = 3;
    int nbBars = 2;
    double a = 2;
    double b = 1;
    double E = 70000000000;
    double A = 0.01;
    // Fields preparation
    nodes.resize(nbNodes*2);
    lockings.resize(nbNodes*2);
    bars.resize(nbBars);
    incidence.resize(nbNodes);
    for(int i = 0 ; i<nbNodes ; i++){
        std::vector<int> nodeIncidence;
        incidence[i] = nodeIncidence;
    }
    // Nodal positions
    nodes[0] = 0.0; // x component of node 0
    nodes[1] = 0.0; // y component of node 0
    nodes[2] = a;
    nodes[3] = -b; // TO HAVE POSITIVE VALUES !!!!
    nodes[4] = 2*a;
    nodes[5] = 0.0;
    // Bars definition
    bars[0] = createBar(0, 1, nodes[0], nodes[1], nodes[2], nodes[3], E, A);
    bars[1] = createBar(2, 1, nodes[4], nodes[5], nodes[2], nodes[3], E, A);
    // Incidence description
    incidence[0].push_back(0);
    incidence[1].push_back(0);
    incidence[1].push_back(1);
    incidence[2].push_back(1);
    // Constraints
    lockings[0] = 1;
    lockings[1] = 1;
    lockings[2] = 1;
    lockings[3] = 0;
    lockings[4] = 1;
    lockings[5] = 1;
    // Number of DOF determination
    int nbDof = 0;
    for(int i=0 ; i<lockings.size() ; i++){
        if(!lockings[i])
            nbDof++;
    }
    // Truss filling
    truss.nodes = nodes;
    truss.incidence = incidence;
    truss.bars = bars;
    truss.lockings = lockings;
    truss.nbDof = nbDof;
    // Return
    return truss;
}

Truss threeBarTruss(){
    Truss truss;
    // Field declaration
    std::vector<double> nodes;
    std::vector<Bar> bars;
    std::vector<std::vector<int> > incidence;
    std::vector<bool> lockings;
    // Field declaration
    int nbNodes = 4;
    int nbBars = 3;
    // Parameters
    double a = 2;
    double b = 1;
    double E1 = 70000000000;
    double E2 = 5000000000;
    double A = 0.01;
    // Fields preparation
    nodes.resize(nbNodes*2);
    lockings.resize(nbNodes*2);
    bars.resize(nbBars);
    incidence.resize(nbNodes);
    for(int i = 0 ; i<nbNodes ; i++){
        std::vector<int> nodeIncidence;
        incidence[i] = nodeIncidence;
    }
    // Nodal positions
    nodes[0] = 0.0; // x component of node 0
    nodes[1] = 0.0; // y component of node 0
    nodes[2] = a;
    nodes[3] = -b;
    nodes[4] = 2*a;
    nodes[5] = 0.0;
    nodes[6] = a;
    nodes[7] = -b - sqrt(a*a+b*b);
    // Bars definition
    bars[0] = createBar(0, 1, nodes[0], nodes[1], nodes[2], nodes[3], E1, A);
    bars[1] = createBar(2, 1, nodes[4], nodes[5], nodes[2], nodes[3], E1, A);
    bars[2] = createBar(1, 3, nodes[2], nodes[3], nodes[6], nodes[7], E2, A);
    // Incidence description
    incidence[0].push_back(0);
    incidence[1].push_back(0);
    incidence[1].push_back(1);
    incidence[1].push_back(2);
    incidence[2].push_back(1);
    incidence[3].push_back(2);
    // Constraints
    lockings[0] = 1;
    lockings[1] = 1;
    lockings[2] = 1;
    lockings[3] = 0;
    lockings[4] = 1;
    lockings[5] = 1;
    lockings[6] = 1;
    lockings[7] = 0;
    // Number of DOF determination
    int nbDof = 0;
    for(int i=0 ; i<lockings.size() ; i++){
        if(!lockings[i])
            nbDof++;
    }
    // Truss filling
    truss.nodes = nodes;
    truss.incidence = incidence;
    truss.bars = bars;
    truss.lockings = lockings;
    truss.nbDof = nbDof;
    // Return
    return truss;
}

Bar createBar(int nodeOne, int nodeTwo, double nodeOneX, double nodeOneY,
    double nodeTwoX, double nodeTwoY, double E, double A){
    // Bar declaration
    Bar bar;
    // Field initialization
    bar.E = E;
    bar.A = A;
    bar.nodes[0] = nodeOne;
    bar.nodes[1] = nodeTwo;
    bar.l0 = sqrt(power(nodeOneX-nodeTwoX, 2) + power(nodeOneY-nodeTwoY, 2));
    // Return
    return bar;
}

// For a simple 1 DOF case
double PVW(Truss &truss, double u, double F){
    std::vector<double> uVec(1,u);
    std::vector<double> FVec(1,F);
    std::vector<double> OOBVec(1);
    PVW(truss, uVec, FVec, OOBVec);
    return OOBVec[0];
}

void PVW(Truss &truss, std::vector<double> &u, std::vector<double> &F, double lambda,
    std::vector<double> &OOB){
    // Call the usual function
    std::vector<double> realForce(truss.nbDof);
    sv(lambda, F, realForce);
    PVW(truss, u, realForce, OOB);
}

void PVW(Truss &truss, std::vector<double> &u, std::vector<double> &F,
    std::vector<double> &OOB){
    // Creates a complete uC vector from u with value 0 if not a DOF
    std::vector<double> uC(truss.nodes.size());
    int dof = 0;
    for(int i=0 ; i<truss.nodes.size() ; i++){
        if(truss.lockings[i])
            uC[i] = 0.0;
        else{
            uC[i] = u[dof];
            dof++;
        }
    }

    // Calculation of dWInt and dWExt
    dof = 0;
    for(int i=0 ; i < truss.nodes.size() ; i++){
        if(!truss.lockings[i]){
            OOB[dof] = 0.0;
            // dWInt
            int node = i/2; // integer division
            // Loop only over the bars that are incident to the node
            // (the others are not affected by the virtual displacement)
            for(int barNb=0 ; barNb<truss.incidence[node].size() ; barNb++){
                Bar bar = truss.bars[truss.incidence[node][barNb]];
                // a and b determination (with sign)
                double a,b;
                int otherNode;
                if(bar.nodes[1] == node)
                    otherNode = 0;
                else
                    otherNode = 1;
                a = truss.nodes[2*node] - truss.nodes[2*bar.nodes[otherNode]];
                b = truss.nodes[2*node+1] - truss.nodes[2*bar.nodes[otherNode]+1];
                // Current length
                double l2;
                if(!(i%2)){ // DOF along x-direction
                        l2 = power(a+uC[i]-uC[2*bar.nodes[otherNode]],2) + power(b+uC[i+1]-uC[2*bar.nodes[otherNode]+1],2);
                }
                else{ // DOF along y-direction
                        l2 = power(a+uC[i-1]-uC[2*bar.nodes[otherNode]],2) + power(b+uC[i]-uC[2*bar.nodes[otherNode]+1],2);
                }
                // Green-Lagrange virtual strain
                double dE11;
                if(i%2){dE11 = (b+uC[i]-uC[2*bar.nodes[otherNode]+1])/power(bar.l0,2);}
                else{dE11 = (a+uC[i]-uC[2*bar.nodes[otherNode]])/power(bar.l0,2);}
                // Piola stress (for the given consitutive law !!!)
                double S11 = bar.E/2.0 * (l2 - power(bar.l0,2))/power(bar.l0,2);
                // dWInt update
                OOB[dof] += bar.A * bar.l0 * S11 * dE11;
            }
            // OOB = dWInt - dWExt
            OOB[dof] -= F[dof];
            // Increment DOF
            dof++;
        }
    }
}
