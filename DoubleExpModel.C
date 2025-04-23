#include <iostream>
#include <fstream>
#include <TMath.h>


Double_t deFunc(Double_t *x, Double_t *p){
    return p[0]*TMath::Exp((-1)*x[0]/p[1])+p[2]*TMath::Exp((-1)*x[0]/p[3]);
}

/*
PARAMETERS of the double exponent model:
p[0] is the constant term for \lambda_1
p[1] -> \lambda_1
p[2] is the constant term for \lambda_2
p[3] -> \lambda_2
*/

Double_t deModel(std::string filename = "250410"){
    //opening the file with fitted spectra
    // TFile *spectra = new TFile(filename.c_str(), "READ");
    // if(!spectra->IsOpen()){
    //     std::cerr << "File not opened!" << std::endl;
    //     return -1;
    // }

    //opening the output file
    TFile *output = new TFile(Form("%s.root", filename.c_str()), "RECREATE");
    // if(!spectra->IsOpen()){
    //     std::cerr << "File not opened!" << std::endl;
    //     return -1;
    // }

    //creating the graph and the function for fitting
    TGraphErrors *DE = new TGraphErrors(Form("%s.txt", filename.c_str()), "%lg %lg %lg %lg");

    TF1 *deF = new TF1("deF", deFunc, 0, 300, 4);
    deF->SetParNames("I1", "lambda1", "I2", "lambda2");
    deF->SetParameters(3, 450, 3, 40);
    deF->SetParLimits(0, 0, 1e2);
    deF->SetParLimits(1, 0, 1e3);
    deF->SetParLimits(2, 0, 1e2);
    deF->SetParLimits(3, 0, 400);
   
    //reading data from the file with fitted spectra
    // for(){
    //     //dostać się do Qh;
    // }

    //fitting the graph
    TFitResultPtr fit = DE->Fit(deF, "S");

    std::string paramsString = Form("#splitline{#splitline{I_{1} = (%.3f +/- %.3f) a.u.}{#lambda_{1} = (%.2f +/- %.2f) cm}}{#splitline{I_{2} = (%.2f +/- %.2f) a.u.}{#lambda_{2} = (%.2f +/- %.2f) cm}}", deF->GetParameter(0), deF->GetParError(0), 
                                            deF->GetParameter(1), deF->GetParError(1), deF->GetParameter(2), deF->GetParError(2), deF->GetParameter(3), deF->GetParError(3));
    TLatex *paramsText = new TLatex();

    //creating the canvas for the graph
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    DE->SetTitle("DE Model; position [cm]; Q [C]");
    DE->SetMarkerStyle(20);
    DE->SetMarkerSize(0.5);
    DE->SetMarkerColor(kBlue);
    DE->SetLineColor(kBlue);
    DE->SetLineWidth(2);
    DE->Draw("AP");
    paramsText->DrawLatex(.8, .8, paramsString.c_str());

    c1->Write();
    DE->Write();

    return 0;
}