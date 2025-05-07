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

Double_t singFunc(Double_t *x, Double_t *p){
    return p[0]*TMath::Exp((-1)*x[0]/p[1])+p[2];
}

Double_t Model(std::string directory = "/home/szymon/LHCb/20250428firstMeasCitiroc/BCF20XL1/"){
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

    std::vector <Float_t> mean, erMean, pos, erPos;
    for(Int_t i = 1; i<41; i++){
        std::string hName = Form("Run%i_PHA_LG_0_50.txt", i);

        TH1I *hist = (TH1I*)spectra->Get(hName.c_str());
        std::string fName = Form("%s_fit", hName.c_str());
        TF1 *func = hist->GetFunction(fName.c_str());

        Float_t meanS = func->GetParameter(1);
        Float_t erMeanS = func->GetParError(1);
        std::cout<<meanS<<std::endl;
        mean.push_back(meanS);
        erMean.push_back(erMeanS);
        pos.push_back(i*10);
        erPos.push_back(0.1);
    }
    Int_t n = pos.size();

    //creating the graphs and the functions for fitting
    TGraphErrors *DE = new TGraphErrors(n, pos.data(), mean.data(), erPos.data(), erMean.data());
    TGraphErrors *SE = new TGraphErrors(*DE);

    TF1 *deF = new TF1("deF", deFunc, 0, 410, 4);
    deF->SetParNames("I1", "lambda1", "I2", "lambda2");
    deF->SetParameters(3, 4500, 3, 40);
    deF->SetParLimits(0, 0, 1e2);
    //deF->SetParLimits(1, 0, 6e3);
    deF->SetParLimits(2, 0, 1e2);
    deF->SetParLimits(3, 0, 40);

    TF1 *siF = new TF1("siF", singFunc, 0, 410, 3);
    siF->SetParNames("I0", "lambda", "const");
    siF->SetParameters(30, 450, 30);

   
    //reading data from the file with fitted spectra
    // for(){
    //     //dostać się do Qh;
    // }

    //fitting the graphs
    TFitResultPtr fit = DE->Fit(deF, "S");
    TFitResultPtr fitSE = SE->Fit(siF, "S");

    // Component 1: [0]*exp(-x/[1])
    TF1 *f1 = new TF1("f1", "[0]*exp(-x/[1])", 0, 410);
    f1->SetParameters(deF->GetParameter(0), deF->GetParameter(1));
    f1->SetLineColor(kOcean);

    // Component 2: [2]*exp(-x/[3])
    TF1 *f2 = new TF1("f2", "[0]*exp(-x/[1])", 0, 410);
    f2->SetParameters(deF->GetParameter(2), deF->GetParameter(3));
    f2->SetLineColor(kMagenta);

    std::string paramsString = Form("#splitline{#splitline{I_{1} = (%.3f +/- %.3f) a.u.}{#lambda_{1} = (%.2f +/- %.2f) cm}}{#splitline{I_{2} = (%.2f +/- %.2f) a.u.}{#lambda_{2} = (%.2f +/- %.2f) cm}}", deF->GetParameter(0), deF->GetParError(0), 
                                            deF->GetParameter(1), deF->GetParError(1), deF->GetParameter(2), deF->GetParError(2), deF->GetParameter(3), deF->GetParError(3));
    TLatex *paramsText = new TLatex();

    TLegend *legDE = new TLegend();

    //creating the canvas for the graph
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 600);
    c1->Divide(2, 1);

    c1->cd(1);
    DE->SetTitle("DE Model; position [cm]; Q [a.u.]");
    DE->SetMarkerStyle(20);
    DE->SetMarkerSize(0.5);
    DE->SetMarkerColor(kBlue);
    DE->SetLineColor(kBlue);
    DE->SetLineWidth(2);
    DE->Draw("AP");

    f1->SetLineStyle(2);
    f2->SetLineStyle(2);
    legDE->AddEntry(f1, "long component", "l");
    legDE->AddEntry(f2, "short component", "l");
    legDE->AddEntry(deF, "DE model", "l");

    f1->Draw("SAME");
    f2->Draw("SAME");
    paramsText->DrawLatex(300, 35, paramsString.c_str());
    legDE->Draw();

    c1->cd(2);
    SE->SetTitle("SE Model; position [cm]; Q [a.u.]");
    SE->SetMarkerStyle(20);
    SE->SetMarkerSize(0.5);
    SE->SetMarkerColor(kBlue);
    SE->SetLineColor(kBlue);
    SE->SetLineWidth(2);
    SE->Draw("AP");

    c1->Write();
    DE->Write();

    return 0;
}