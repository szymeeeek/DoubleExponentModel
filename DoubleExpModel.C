#include <iostream>
#include <fstream>
#include <TMath.h>


Double_t deFunc(Double_t *x, Double_t *p){
    return p[0]*TMath::Exp((-1)*x[0]/p[1])+p[2]*TMath::Exp((-1)*x[0]/p[3])+p[4];
}

/*
PARAMETERS of the double exponent model:
p[0] is the constant term for \lambda_1
p[1] -> \lambda_1
p[2] is the constant term for \lambda_2
p[3] -> \lambda_2
p[4] - constant
*/

Double_t singFunc(Double_t *x, Double_t *p){
    return p[0]*TMath::Exp((-1)*x[0]/p[1])+p[2];
}

Double_t ModelCitiroc(std::string directory = "/home/szymon/LHCb/20250524testsRep/BCF20XL1/Citiroc/"){
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
    pos = {3, 20, 40, 60, 80, 100, 150, 200, 300, 400};

    for(Int_t i = 1; i<42; i++){
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

    TF1 *deF = new TF1("deF", deFunc, 0, 410, 5);
    deF->SetParNames("I1", "lambda1", "I2", "lambda2");
    deF->SetParameters(35, 450, 3, 40, 30);
    deF->SetParLimits(0, 0, 1e2);
    deF->SetParLimits(1, 0, 6e3);
    deF->SetParLimits(2, 0, 1e2);
    deF->SetParLimits(3, 0, 100);
    deF->SetParLimits(4, 0, 30);

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

Double_t ModelOsc(std::string directory = "/home/szymon/LHCb/20250524testsRep/BCF20XL1/Scope/"){
    Bool_t debug = kFALSE;
    //opening the output file
    std::string outFilename = Form("%soutput_modelFit.root", directory.c_str());
    TFile *output = new TFile(outFilename.c_str(), "RECREATE");
    if(!output->IsOpen()){
        std::cerr << "File not opened!" << std::endl;
        return -1;
    }

    std::vector <Float_t> mean, erMean, pos, erPos;
    //----------extracting data from files with fitted spectra----------
    // List of input ROOT files
    Int_t firstFile = 20;
    Int_t lastFile = 31;

    pos = {3, 20, 40, 60, 80, 100, 150, 200, 300, 400};
    // Loop over files
    for (Int_t i = firstFile; i <= lastFile; ++i) {
        // Construct file name and histogram name
        std::string fileName = Form("%s_scope_%i_HISTOS.root", directory.c_str(), i);
        TFile* file = TFile::Open(fileName.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Cannot open file: " << fileName << std::endl;
            continue;
        }

        // Get histogram
        TH1* hist = dynamic_cast<TH1*>(file->Get("Qh"));
        if (!hist) {
            std::cerr << "Histogram not found in file: " << fileName << std::endl;
            file->Close();
            continue;
        }

        // Fit with Gaussian
        TF1 *func = hist->GetFunction("QGaus");

        // Extract fit parameters
        mean.push_back(func->GetParameter(1));    
        erMean.push_back(func->GetParError(1));   
        erPos.push_back(0.1); // Assuming constant error for position
        if(pos.empty()) {
            pos.push_back(i * 10); // Assuming position is 10*i for each file
        }

        if(debug){
            std::cout << "File: " << fileName << std::endl;
            std::cout << "  Mean  = " << mean.back() << std::endl;
            std::cout << "  MeanEr = " << erMean.back() << std::endl;
        }

        // Clean up
        delete func;
        file->Close();
    }

    Int_t n = pos.size();

    //creating the graphs and the functions for fitting
    TGraphErrors *DE = new TGraphErrors(n, pos.data(), mean.data(), erPos.data(), erMean.data());
    TGraphErrors *SE = new TGraphErrors(*DE);

    TF1 *deF = new TF1("deF", deFunc, 0, 410, 5);
    deF->SetParNames("I1", "lambda1", "I2", "lambda2");
    deF->SetParameters(35, 450, 3, 40, 30);
    deF->SetParLimits(0, 0, 1e2);
    deF->SetParLimits(1, 0, 6e3);
    deF->SetParLimits(2, 0, 1e2);
    deF->SetParLimits(3, 0, 100);
    deF->SetParLimits(4, 0, 30);

    TF1 *siF = new TF1("siF", singFunc, 0, 410, 3);
    siF->SetParNames("I0", "lambda", "const");
    siF->SetParameters(30, 450, 30);

    //fitting the graphs
    TFitResultPtr fit = DE->Fit(deF, "S");
    TFitResultPtr fitSE = SE->Fit(siF, "S");

    //---------components for the single exponential model---------
    // Component 1: [0]*exp(-x/[1])
    TF1 *f1SI = new TF1("f1", "[0]*exp(-x/[1])", 0, 410);
    f1SI->SetParameters(siF->GetParameter(0), siF->GetParameter(1));
    f1SI->SetLineColor(kOcean);
    // Component 2: constant
    TF1 *f2SI = new TF1("f2", "[0]", 0, 410);
    f2SI->SetParameters(siF->GetParameter(2));
    f2SI->SetLineColor(kMagenta);

    //---------components for the double exponential model---------
    // Component 1: [0]*exp(-x/[1])
    TF1 *f1 = new TF1("f1", "[0]*exp(-x/[1])", 0, 410);
    f1->SetParameters(deF->GetParameter(0), deF->GetParameter(1));
    f1->SetLineColor(kOcean);

    // Component 2: [2]*exp(-x/[3])
    TF1 *f2 = new TF1("f2", "[0]*exp(-x/[1])", 0, 410);
    f2->SetParameters(deF->GetParameter(2), deF->GetParameter(3));
    f2->SetLineColor(kMagenta);

    // Component 3: constant
    TF1 *f3 = new TF1("f3", "[0]", 0, 410);
    f3->SetParameter(0, deF->GetParameter(4));
    f3->SetLineColor(kRed);
    f3->SetLineStyle(2);
    //--------------------------------------------------------------

    std::string paramsStringDE = Form("#splitline{#splitline{I_{1} = (%.3f +/- %.3f) a.u.}{#lambda_{1} = (%.2f +/- %.2f) cm}}{#splitline{I_{2} = (%.2f +/- %.2f) a.u.}{#lambda_{2} = (%.2f +/- %.2f) cm}}", deF->GetParameter(0), deF->GetParError(0), 
                                            deF->GetParameter(1), deF->GetParError(1), deF->GetParameter(2), deF->GetParError(2), deF->GetParameter(3), deF->GetParError(3));
    std::string paramsStringSE = Form("#splitline{#splitline{I_{0} = (%.3f +/- %.3f) a.u.}{#lambda = (%.2f +/- %.2f) cm}}{#splitline{const = (%.2f +/- %.2f) a.u.}}", siF->GetParameter(0), siF->GetParError(0), 
                                            siF->GetParameter(1), siF->GetParError(1), siF->GetParameter(2), siF->GetParError(2));
    
    
    TLatex *paramsTextDE = new TLatex();
    TLatex *paramsTextSE = new TLatex();

    TLegend *legDE = new TLegend();
    TLegend *legSE = new TLegend();

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
    legDE->AddEntry(f3, "constant", "l");
    legDE->AddEntry(deF, "DE model", "l");

    f1->Draw("SAME");
    f2->Draw("SAME");
    paramsTextDE->DrawLatex(300, 35, paramsStringDE.c_str());
    legDE->Draw();

    c1->cd(2);
    SE->SetTitle("SE Model; position [cm]; Q [a.u.]");
    SE->SetMarkerStyle(20);
    SE->SetMarkerSize(0.5);
    SE->SetMarkerColor(kBlue);
    SE->SetLineColor(kBlue);
    SE->SetLineWidth(2);
    SE->Draw("AP");

    f1->SetLineStyle(2);
    f2->SetLineStyle(2);
    legSE->AddEntry(f1, "single component", "l");
    legSE->AddEntry(f2, "constant", "l");
    legSE->AddEntry(siF, "SE model", "l");
    f1->Draw("SAME");
    f2->Draw("SAME");
    paramsTextSE->DrawLatex(300, 35, paramsStringSE.c_str());
    legSE->Draw();

    c1->Write();
    DE->Write();

    return 0;
}