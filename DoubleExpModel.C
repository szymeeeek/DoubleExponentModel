#include <iostream>
#include <fstream>
#include <TMath.h>

//changelog
// 2024-06-10: Added the chi2/ndf calculation for both models.

//This macro provides the functions for fitting the double and sinle exponential models to the data from Citiroc and Oscilloscope.
//To use the following macro, you need to first process the data with the Citiroc or Oscilloscope analysis macros.

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

/*
PARAMETERS of the double exponent model:
p[0] is the constant term for \lambda
p[1] -> \lambda
p[2] - constant
*/

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
    
    //if the measurement was performed with a fixed step, leave the pos vector empty
    pos = {3, 20, 40, 60, 80, 100, 150, 200, 300, 400};
    Double_t step = 10; //cm

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
        if(pos.empty()){
            pos.push_back(i*step);
        }
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
    siF->SetParameters(16, 45, 0);
    siF->SetParLimits(0, 0, 1e2);
    siF->SetParLimits(1, 0, 6e2);
    siF->SetParLimits(2, 0, 14);

    //fitting the graphs
    TFitResultPtr fit = DE->Fit(deF, "R");
    TFitResultPtr fitSE = SE->Fit(siF, "R");

    Double_t chi2DE = deF->GetChisquare();
    Int_t ndfDE = deF->GetNDF();
    Double_t chi2PerNdfDE = (ndfDE != 0) ? chi2DE / ndfDE : 0;

    Double_t chi2SE = siF->GetChisquare();
    Int_t ndfSE = siF->GetNDF();
    Double_t chi2PerNdfSE = (ndfSE != 0) ? chi2SE / ndfSE : 0;

    std::cout << "Double Exp: chi2 = " << chi2DE << ", ndf = " << ndfDE << ", chi2/ndf = " << chi2PerNdfDE << std::endl;
    std::cout << "Single Exp: chi2 = " << chi2SE << ", ndf = " << ndfSE << ", chi2/ndf = " << chi2PerNdfSE << std::endl;

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

    std::string paramsStringDE = Form(
        "#splitline{#splitline{I_{1} = (%.2f #pm %.2f) a.u.}{#lambda_{1} = (%.1f #pm %.1f) cm}}"
        "{#splitline{I_{2} = (%.1f #pm %.1f) a.u.}{#splitline{#lambda_{2} = (%.2f #pm %.2f) cm}{const. = (%.4f #pm %.4f) a.u.}}}"
        "\\n#chi^{2}/ndf = %.2f",
        deF->GetParameter(0), deF->GetParError(0), 
        deF->GetParameter(1), deF->GetParError(1), 
        deF->GetParameter(2), deF->GetParError(2), 
        deF->GetParameter(3), deF->GetParError(3), 
        deF->GetParameter(4), deF->GetParError(4),
        chi2PerNdfDE
    );
    std::string paramsStringSE = Form(
        "#splitline{#splitline{I_{0} = (%.2f #pm %.2f) a.u.}{#lambda = (%.2f #pm %.2f) cm}}"
        "{#splitline{const. = (%.4f #pm %.4f) a.u.}{}}"
        "\\n#chi^{2}/ndf = %.2f",
        siF->GetParameter(0), siF->GetParError(0), 
        siF->GetParameter(1), siF->GetParError(1), 
        siF->GetParameter(2), siF->GetParError(2),
        chi2PerNdfSE
    );

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
    DE->SetMinimum(0);
    DE->SetMaximum(50);
    DE->Draw("AP");

    f1->SetLineStyle(2);
    f2->SetLineStyle(2);
    legDE->AddEntry(f1, "long component", "l");
    legDE->AddEntry(f2, "short component", "l");
    legDE->AddEntry(f3, "constant", "l");
    legDE->AddEntry(deF, "DE model", "l");

    f1->Draw("SAME");
    f2->Draw("SAME");
    f3->Draw("SAME");
    paramsTextDE->SetNDC(kTRUE);
    paramsTextDE->SetTextSize(0.04);
    paramsTextSE->SetTextSize(0.04);
    paramsTextDE->DrawLatex(0.5, 0.75, paramsStringDE.c_str());
    legDE->Draw();

    c1->cd(2);
    SE->SetTitle("SE Model; position [cm]; Q [a.u.]");
    SE->SetMarkerStyle(20);
    SE->SetMarkerSize(0.5);
    SE->SetMarkerColor(kBlue);
    SE->SetLineColor(kBlue);
    SE->SetLineWidth(2);
    SE->SetMinimum(0);
    SE->SetMaximum(50);
    SE->Draw("AP");

    f1SI->SetLineStyle(2);
    f2SI->SetLineStyle(2);
    legSE->AddEntry(f1, "single component", "l");
    legSE->AddEntry(f2, "constant", "l");
    legSE->AddEntry(siF, "SE model", "l");
    f1SI->Draw("SAME");
    f2SI->Draw("SAME");
    paramsTextSE->SetNDC(kTRUE);
    paramsTextSE->DrawLatex(0.5, 0.75, paramsStringSE.c_str());
    legSE->Draw();

    c1->Update();
    output->cd();
    c1->Write();
    DE->Write();
    SE->Write();

    return 0;
}

Double_t ModelOsc(std::string directory = "/scratch3/lhcb/data/20250601testsWithScopeRep/BCF12XL/20250602"){
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
    Int_t firstFile = 22;
    Int_t lastFile = 33;

    pos = {10, 20, 35, 50, 75, 100, 150, 200, 250, 300, 350, 400}; // Assuming positions in cm
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
        TH1D* hist = dynamic_cast<TH1D*>(file->Get("Qh"));
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
    siF->SetParameters(16, 45, 0);
    siF->SetParLimits(0, 0, 1e2);
    siF->SetParLimits(1, 0, 6e2);
    siF->SetParLimits(2, 0, 14);

    //fitting the graphs
    TFitResultPtr fit = DE->Fit(deF, "R");
    TFitResultPtr fitSE = SE->Fit(siF, "R");

    Double_t chi2DE = deF->GetChisquare();
    Int_t ndfDE = deF->GetNDF();
    Double_t chi2PerNdfDE = (ndfDE != 0) ? chi2DE / ndfDE : 0;

    Double_t chi2SE = siF->GetChisquare();
    Int_t ndfSE = siF->GetNDF();
    Double_t chi2PerNdfSE = (ndfSE != 0) ? chi2SE / ndfSE : 0;

    std::cout << "Double Exp: chi2 = " << chi2DE << ", ndf = " << ndfDE << ", chi2/ndf = " << chi2PerNdfDE << std::endl;
    std::cout << "Single Exp: chi2 = " << chi2SE << ", ndf = " << ndfSE << ", chi2/ndf = " << chi2PerNdfSE << std::endl;

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

    std::string paramsStringDE = Form(
        "#splitline"
        "{#splitline{I_{1} = (%.2f #pm %.2f) a.u.}{#lambda_{1} = (%.1f #pm %.1f) cm}}"
        "{#splitline{I_{2} = (%.1f #pm %.1f) a.u.}{#splitline{#lambda_{2} = (%.2f #pm %.2f) cm}{#splitline{const. = (%.4f #pm %.4f) a.u.}{#chi^{2}/ndf = %.2f}}}}",
        deF->GetParameter(0), deF->GetParError(0), 
        deF->GetParameter(1), deF->GetParError(1), 
        deF->GetParameter(2), deF->GetParError(2), 
        deF->GetParameter(3), deF->GetParError(3), 
        deF->GetParameter(4), deF->GetParError(4),
        chi2PerNdfDE
    );
    std::string paramsStringSE = Form(
        "#splitline{#splitline{I_{0} = (%.2f #pm %.2f) a.u.}{#lambda = (%.2f #pm %.2f) cm}}"
        "{#splitline{const. = (%.4f #pm %.4f) a.u.}{#chi^{2}/ndf = %.2f}}",
        siF->GetParameter(0), siF->GetParError(0), 
        siF->GetParameter(1), siF->GetParError(1), 
        siF->GetParameter(2), siF->GetParError(2),
        chi2PerNdfSE
    );
    
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
    DE->SetMinimum(0);
    DE->SetMaximum(50);
    DE->Draw("AP");

    f1->SetLineStyle(2);
    f2->SetLineStyle(2);
    legDE->AddEntry(f1, "long component", "l");
    legDE->AddEntry(f2, "short component", "l");
    legDE->AddEntry(f3, "constant", "l");
    legDE->AddEntry(deF, "DE model", "l");

    f1->Draw("SAME");
    f2->Draw("SAME");
    f3->Draw("SAME");
    paramsTextDE->SetNDC(kTRUE);
    paramsTextDE->SetTextSize(0.04);
    paramsTextSE->SetTextSize(0.04);
    paramsTextDE->DrawLatex(0.5, 0.75, paramsStringDE.c_str());
    legDE->Draw();

    c1->cd(2);
    SE->SetTitle("SE Model; position [cm]; Q [a.u.]");
    SE->SetMarkerStyle(20);
    SE->SetMarkerSize(0.5);
    SE->SetMarkerColor(kBlue);
    SE->SetLineColor(kBlue);
    SE->SetLineWidth(2);
    SE->SetMinimum(0);
    SE->SetMaximum(50);
    SE->Draw("AP");

    f1SI->SetLineStyle(2);
    f2SI->SetLineStyle(2);
    legSE->AddEntry(f1, "single component", "l");
    legSE->AddEntry(f2, "constant", "l");
    legSE->AddEntry(siF, "SE model", "l");
    f1SI->Draw("SAME");
    f2SI->Draw("SAME");
    paramsTextSE->SetNDC(kTRUE);
    paramsTextSE->DrawLatex(0.5, 0.75, paramsStringSE.c_str());
    legSE->Draw();

    c1->Update();
    output->cd();
    c1->Write();
    DE->Write();
    SE->Write();

    return 0;
}