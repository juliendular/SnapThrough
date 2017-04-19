#include "main.h"

/* Computes the integer power of a given base */
double power(double base, int exponent){
    double result = 1.0;
    for(int i=0 ; i<exponent ; i++){
        result = result * base;
    }
    return result;
}
/* Computes the sum of two vectors */
void vvPlus(std::vector<double> &v1, std::vector<double> &v2, std::vector<double> &sum){
    sum.resize(v1.size());
    for(int i=0 ; i<v1.size() ; i++)
        sum[i] = v1[i] + v2[i];
    return;
}
/* Computes a scalar-vector product */
void sv(double s, std::vector<double> &v, std::vector<double> &p){
    p.resize(v.size());
    for(int i=0 ; i<v.size() ; i++)
        p[i] = s*v[i];
    return;
}
/* Computes a vector-vector product */
double vv(std::vector<double> &v1, std::vector<double> &v2){
    double scalar = 0.0;
    for(int i=0 ; i<v1.size() ; i++)
        scalar += v1[i]*v2[i];
    return scalar;
}
/* Computes a matrix-vector product */
void mv(std::vector<std::vector<double> > &m, std::vector<double> &v,
    std::vector<double> &p){
    p.resize(m.size());
    for(int i=0 ; i<m.size() ; i++){
        double element = 0.0;
        for(int j=0 ; j<v.size() ; j++)
            element +=m[i][j]*v[j];
        p[i] = element;
    }
}
/* Computes the inverse of a square matrix */
int inv(std::vector<std::vector<double> > &m, std::vector<std::vector<double> > &invM){
    if(m.size()==1){
        if(m[0][0] == 0.0){
            std::cout << "Division by zero" << std::endl;
            return 1;
        }
        else{
            invM[0][0] = 1.0/m[0][0];
            return 0;
        }
    }
    else if(m.size()==2){
        double det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        if(det==0){
            std::cout << "Singular matrix" << std::endl;
            return 1;
        }
        else{
            double invDet = 1.0 / det;
            invM[0][0] = invDet * m[1][1];
            invM[0][1] = - invDet * m[0][1];
            invM[1][0] = - invDet * m[1][0];
            invM[1][1] = invDet * m[0][0];
            return 0;
        }
    }
    else{
        std::cout << "n*n (n>2) matrix inversion not implemented" << std::endl;
        return 1;
    }
}
/* Determines if a matrix if positive definite. Actually gives the sign of its determinant */
bool posDefinite(std::vector<std::vector<double> > &m){
    if(m.size()==1){
        return (m[0][0]>0); // Strictly
    }
    else if(m.size()==2){
        //double trace = m[0][0] + m[1][1];
        double det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        //double eig1 = (trace + sqrt(power(trace,2)-4*det))/2.0;
        //double eig2 = (trace - sqrt(power(trace,2)-4*det))/2.0;
        //return (eig1>0 && eig2>0); // BETTER TO USE THE DETERMINANT (see references [2] and [3])
        return det>0; // The name of the function is not appropriate
    }
    else{
        std::cout << "n*n (n>2) positive definitess function not implemented" << std::endl;
        return false;
    }
}
