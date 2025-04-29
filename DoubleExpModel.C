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

Double_t deModel(std::string directory = "/home/szymon/LHCb/20250428firstMeasCitiroc/BCF20XL2/"){
    //opening the file with fitted spectra
    std::string specFile = Form("%soutput_HISTOS.root", directory.c_str());
    TFile *spectra = new TFile(specFile.c_str(), "READ");
    if(!spectra->IsOpen()){
        std::cerr << "File not opened!" << std::endl;
        return -1;
    }

    //opening the output file
    TFile *output = new TFile("output_DEmodel.root", "RECREATE");
    if(!output->IsOpen()){
        std::cerr << "File not opened!" << std::endl;
        return -1;
    }

    Float_t m = 9.0/35.0;
    std::cout<<m<<std::endl;
    std::vector <Float_t> mean, erMean, pos, erPos;
    for(Int_t i = 1; i<41; i++){
        std::string hName = Form("Run%i_PHA_LG_0_50.txt", i);

        TH1I *hist = (TH1I*)spectra->Get(hName.c_str());
        std::string fName = Form("%s_fit", hName.c_str());
        TF1 *func = hist->GetFunction(fName.c_str());

        Float_t meanS = func->GetParameter(1);
        Float_t erMeanS = func->GetParError(1);
        std::cout<<meanS<<std::endl;
        mean.push_back(m*meanS);
        erMean.push_back(m*erMeanS);
        std::cout<<mean.at(i-1)<<" "<<erMean.at(i-1)<<std::endl;
        pos.push_back(i*10);
        erPos.push_back(0.1);
    }
    Int_t n = pos.size();

    //creating the graph and the function for fitting
    //TGraphErrors *DE = new TGraphErrors(Form("%s.txt", filename.c_str()), "%lg %lg %lg %lg");
    TGraphErrors *DE = new TGraphErrors(n, pos.data(), mean.data(), erPos.data(), erMean.data());

    TF1 *deF = new TF1("deF", deFunc, 50, 410, 4);
    deF->SetParNames("I1", "lambda1", "I2", "lambda2");
    deF->SetParameters(3, 4500, 3, 40);
    deF->SetParLimits(0, 0, 1e2);
    //deF->SetParLimits(1, 0, 6e3);
    deF->SetParLimits(2, 0, 1e2);
    deF->SetParLimits(3, 0, 40);
   
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
    DE->SetTitle("DE Model; position [cm]; Q [a.u.]");
    DE->SetMarkerStyle(20);
    DE->SetMarkerSize(0.5);
    DE->SetMarkerColor(kBlue);
    DE->SetLineColor(kBlue);
    DE->SetLineWidth(2);
    DE->Draw("AP");
    paramsText->DrawLatex(300, 35, paramsString.c_str());

    c1->Write();
    DE->Write();

    return 0;
}