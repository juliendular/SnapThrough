#include "main.h"

/* Defines the geometry of the simple truss */
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

/* Defines the geometry of a three-bar truss */
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
    double E2 = 70000000000;
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
    nodes[3] = -b; // To meet the positive conventions of the schematics
    nodes[4] = 2*a;
    nodes[5] = 0.0;
    nodes[6] = a;
    nodes[7] = -b - sqrt(a*a+b*b); // To meet the positive conventions of the schematics
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

Truss fiveBarTruss(){
    Truss truss;
    // Field declaration
    std::vector<double> nodes;
    std::vector<Bar> bars;
    std::vector<std::vector<int> > incidence;
    std::vector<bool> lockings;
    // Field declaration
    int nbNodes = 6;
    int nbBars = 5;
    // Parameters
    double a = 2;
    double b = 1;
    double E1 = 100000000000;
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
    nodes[3] = -b; // To meet the positive conventions of the schematics
    nodes[4] = 2*a;
    nodes[5] = 0.0;

    nodes[6] = 0.0; // x component of node 0
    nodes[7] = 0.0-b; // y component of node 0
    nodes[8] = a;
    nodes[9] = -2*b; // To meet the positive conventions of the schematics
    nodes[10] = 2*a;
    nodes[11] = 0.0-b;
    // Bars definition
    bars[0] = createBar(0, 1, nodes[0], nodes[1], nodes[2], nodes[3], E1, A);
    bars[1] = createBar(2, 1, nodes[4], nodes[5], nodes[2], nodes[3], E1, A);
    bars[2] = createBar(1, 4, nodes[2], nodes[3], nodes[8], nodes[9], E1, A);
    bars[3] = createBar(3, 4, nodes[6], nodes[7], nodes[8], nodes[9], E1, A);
    bars[4] = createBar(4, 5, nodes[8], nodes[9], nodes[10], nodes[11], E1, A);
    // Incidence description
    incidence[0].push_back(0);
    incidence[1].push_back(0);
    incidence[1].push_back(1);
    incidence[1].push_back(2);
    incidence[2].push_back(1);
    incidence[3].push_back(3);
    incidence[4].push_back(2);
    incidence[4].push_back(3);
    incidence[4].push_back(4);
    incidence[5].push_back(4);
    // Constraints
    lockings[0] = 1;
    lockings[1] = 1;
    lockings[2] = 1;
    lockings[3] = 0;
    lockings[4] = 1;
    lockings[5] = 1;
    lockings[6] = 1;
    lockings[7] = 1;
    lockings[8] = 1;
    lockings[9] = 0;
    lockings[10] = 1;
    lockings[11] = 1;
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

Truss sevenBarTruss(){
    Truss truss;
    // Field declaration
    std::vector<double> nodes;
    std::vector<Bar> bars;
    std::vector<std::vector<int> > incidence;
    std::vector<bool> lockings;
    // Field declaration
    int nbNodes = 6;
    int nbBars = 7;
    // Parameters
    double a = 2;
    double b = 1;
    double E1 = 100000000000;
    double E2 = 100000000000;
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
    nodes[3] = -b; // To meet the positive conventions of the schematics
    nodes[4] = 2*a;
    nodes[5] = 0.0;

    nodes[6] = 0.0; // x component of node 0
    nodes[7] = 0.0-b; // y component of node 0
    nodes[8] = a;
    nodes[9] = -2*b; // To meet the positive conventions of the schematics
    nodes[10] = 2*a;
    nodes[11] = 0.0-b;
    // Bars definition
    bars[0] = createBar(0, 1, nodes[0], nodes[1], nodes[2], nodes[3], E1, A);
    bars[1] = createBar(2, 1, nodes[4], nodes[5], nodes[2], nodes[3], E1, A);
    bars[2] = createBar(1, 4, nodes[2], nodes[3], nodes[8], nodes[9], E1, A);
    bars[3] = createBar(3, 4, nodes[6], nodes[7], nodes[8], nodes[9], E1, A);
    bars[4] = createBar(4, 5, nodes[8], nodes[9], nodes[10], nodes[11], E1, A);
    bars[5] = createBar(3, 1, nodes[6], nodes[7], nodes[2], nodes[3], E2, A);
    bars[6] = createBar(5, 1, nodes[10], nodes[11], nodes[2], nodes[3], E2, A);
    // Incidence description
    incidence[0].push_back(0);
    incidence[1].push_back(0);
    incidence[1].push_back(1);
    incidence[1].push_back(2);
    incidence[1].push_back(5);
    incidence[1].push_back(6);
    incidence[2].push_back(1);
    incidence[3].push_back(3);
    incidence[3].push_back(5);
    incidence[4].push_back(2);
    incidence[4].push_back(3);
    incidence[4].push_back(4);
    incidence[5].push_back(4);
    incidence[5].push_back(6);
    // Constraints
    lockings[0] = 1;
    lockings[1] = 1;
    lockings[2] = 1;
    lockings[3] = 0;
    lockings[4] = 1;
    lockings[5] = 1;
    lockings[6] = 1;
    lockings[7] = 1;
    lockings[8] = 1;
    lockings[9] = 0;
    lockings[10] = 1;
    lockings[11] = 1;
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


/* Creates a bar */
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

/* Overload to call the general function in a simple 1 DOF case */
double PVW(Truss &truss, double u, double F){
    std::vector<double> uVec(1,u);
    std::vector<double> FVec(1,F);
    std::vector<double> OOBVec(1);
    PVW(truss, uVec, FVec, OOBVec);
    return OOBVec[0];
}

/* Overload to call the general function with qef (F) and lambda */
void PVW(Truss &truss, std::vector<double> &u, std::vector<double> &F, double lambda,
    std::vector<double> &OOB){
    // Call the usual function
    std::vector<double> realForce(truss.nbDof);
    sv(lambda, F, realForce);
    PVW(truss, u, realForce, OOB);
}

/* GENERIC OUT OF BALANCE FORCE CALCULATION, this function works for any truss */
void PVW(Truss &truss, std::vector<double> &u, std::vector<double> &F,
    std::vector<double> &OOB){
    // Creates a complete uC vector from u with value 0 if not a DOF (simplifies the following)
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
    // Calculation of FInt and FExt (the virtual displacement is put in evidence and does not appear (OK, it is arbitrary))
    dof = 0;
    // Loop over the nodal displacement, do something only if the current one is a DOF
    for(int i=0 ; i < truss.nodes.size() ; i++){
        if(!truss.lockings[i]){ // Is the current displacement a DOF?
            OOB[dof] = 0.0; // Initialization to 0
            // FInt
            int node = i/2; // Integer division to get the node number
            // Loop only over the bars that are incident to the node
            // (the others are not affected by the virtual displacement)
            for(int barNb=0 ; barNb<truss.incidence[node].size() ; barNb++){
                // Gets the considered bar
                Bar bar = truss.bars[truss.incidence[node][barNb]];
                // Lengths a and b determination (with sign)
                double a,b; // a is the horizontal length and b the vertical one (with sign conventions...)
                int otherNode;
                if(bar.nodes[1] == node) // On which side of the bar is the current node? (important for sign conventions)
                    otherNode = 0;
                else
                    otherNode = 1;
                a = truss.nodes[2*node] - truss.nodes[2*bar.nodes[otherNode]];
                b = truss.nodes[2*node+1] - truss.nodes[2*bar.nodes[otherNode]+1];
                // Current (deformed) length
                double l2;
                if(!(i%2)) // DOF along x-direction
                        l2 = power(a+uC[i]-uC[2*bar.nodes[otherNode]],2) + power(b+uC[i+1]-uC[2*bar.nodes[otherNode]+1],2);
                else // DOF along y-direction
                        l2 = power(a+uC[i-1]-uC[2*bar.nodes[otherNode]],2) + power(b+uC[i]-uC[2*bar.nodes[otherNode]+1],2);
                // Green-Lagrange virtual strain (analytical derivation similar to what is done in section 1)
                double dE11;
                if(i%2){dE11 = (b+uC[i]-uC[2*bar.nodes[otherNode]+1])/power(bar.l0,2);}
                else{dE11 = (a+uC[i]-uC[2*bar.nodes[otherNode]])/power(bar.l0,2);}
                // Piola stress (for the given consitutive law !!!)
                double S11 = bar.E/2.0 * (l2 - power(bar.l0,2))/power(bar.l0,2);
                // FInt update
                OOB[dof] += bar.A * bar.l0 * S11 * dE11; // Integration over the reference configuration
            }
            // OOB = FInt - FExt
            OOB[dof] -= F[dof]; // Reminder: the virtual displacement is put in evidence
            // Increment DOF
            dof++;
            // All the contributions from the bars adjacent to the considered DOF have been considered, go to the next DOF
        }
    }
}

/* Specific function for the first considered three-bar truss with another constitutive law (sigma=E epsilon), just used as a comparison */
void PVW_cauchy(Truss &truss, std::vector<double> &u, std::vector<double> &F,
    std::vector<double> &OOB){
    // Intermediate variables
    double a = truss.nodes[2];
    double b = -truss.nodes[3];
    double l0 = sqrt(power(a,2)+power(b,2));
    double l = sqrt(power(a,2)+power(b-u[0],2)); // length of the oblique bars
    // OOB forces (see analytical formula in the report)
    OOB[0] = + truss.bars[2].E*truss.bars[2].A / l0 * (u[0]-u[1])
            + truss.bars[0].E*truss.bars[0].A / (l0 * l) * 2 * (u[0]-b)*(l-l0) - F[0];
    OOB[1] = - truss.bars[2].E*truss.bars[2].A / l0 * (u[0]-u[1]) - F[1];

}
