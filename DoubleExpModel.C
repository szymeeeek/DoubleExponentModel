#include <iostream>
#include <fstream>
#include <TMath>


Double_t deFunc(Double_t *x, Double_t *p){
    return p[0]*TMath::Exp((-1)*x/p[1])+p[2]*TMath::Exp((-1)*x/p[3]);
}

/*
PARAMETERS:
p[0] is the constant term for \lambda_1
p[1] -> \lambda_1
p[2] is the constant term for \lambda_2
p[3] -> \lambda_2
*/

Double_t deModel(std::string filename = ""){
    

    return 0;
}