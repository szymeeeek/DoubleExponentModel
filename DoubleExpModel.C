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
    TFile *spectra = new TFile(filename.c_str(), "READ");
    if(!spectra->IsOpen()){
        std::cerr << "File not opened!" << std::endl;
        return -1;
    }

    TGraph *DE = new TGraph();
    TF1 *deF = new TF1("deF", deFunc, 0, 1000, 4);
    deF->SetParNames("I1", "lambda1", "I2", "lambda2");
    deF->SetParameters(1, 1, 1, 1);
    deF->SetParLimits(0, 0, 1e6);
    deF->SetParLimits(1, 0, 1e6);
    deF->SetParLimits(2, 0, 1e6);
    deF->SetParLimits(3, 0, 1e6);
   
    for(){
        //dostać się do Qh;
    }

    TFitResultPtr fit = DE->Fit(deF, "S");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    DE->SetTitle("DE Model; position [mm]; Q [a.u.]");
    DE->SetMarkerStyle(20);
    DE->SetMarkerSize(0.5);
    DE->SetMarkerColor(kBlue);
    DE->SetLineColor(kBlue);
    DE->SetLineWidth(2);

    return 0;
}