#include "main.h"

/* Computes the integer power of a given base
IN: - base: the number to be raised to power
    - exponent: integer power
OUT: the result of base^exponent
*/
double power(double base, int exponent){
    double result = 1.0;
    for(int i=0 ; i<exponent ; i++){
        result = result * base;
    }
    return result;
}
