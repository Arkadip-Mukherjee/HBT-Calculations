
// TO RUN THIS PROGRAM, USE :            g++ -o run correlfit.cpp `root-config --cflags --glibs`

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TH3F.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLegend.h>
#include <TString.h>
#include <TError.h>
#include <Math/MinimizerOptions.h>

using namespace std;

const int q_num_bins = 51;
const double qmin = -0.1;
const double qmax = 0.1;
const int kT_num_bins = 4;
const double kT_min = 0.01;
const double kT_max = 2.1;
double step = (kT_max - kT_min) / kT_num_bins;

// Initial guesses for HBT Radii
const double lambda_ini = 1.0;
const double Rout_ini = 4.0;
const double Rside_ini = 4.0;
const double Rlong_ini = 5.0;
const double Rout_min = 2.0, Rout_max = 10.0;
const double Rside_min = 2.0, Rside_max = 10.0;
const double Rlong_min = 2.0, Rlong_max = 14.0;

Double_t model_func_3D(Double_t *x, Double_t *par) {
	Double_t qoutsq = x[0] * x[0];
	Double_t qsidesq = x[1] * x[1];
	Double_t qlongsq = x[2] * x[2];
	Double_t lambda = TMath::Abs(par[0]);
	Double_t gpart = TMath::Exp((-par[1]*par[1]*qoutsq - par[2]*par[2]*qsidesq - par[3]*par[3]*qlongsq) / 0.038937937);
	return (1 + lambda*gpart);
}


int main() {
	gErrorIgnoreLevel = kWarning;
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
	ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-6);
	
	std::ofstream outputFile("HBT_radii_data_19.6GeV.csv");
	outputFile << "kT,Rout,Rside,Rlong,ParticleType\n"; // Header for the CSV file
	
	std::vector<double> Rout_storeval_1, Rside_storeval_1, Rlong_storeval_1;
	std::vector<double> Rout_storeval_2, Rside_storeval_2, Rlong_storeval_2;
	std::vector<double> Rout_storeval_3, Rside_storeval_3, Rlong_storeval_3;
	std::vector<double> Rout_storeval_4, Rside_storeval_4, Rlong_storeval_4;
	std::vector<double> Rout_storeval_5, Rside_storeval_5, Rlong_storeval_5;
	std::vector<double> Rout_storeval_6, Rside_storeval_6, Rlong_storeval_6;
	
	std::vector<double> ChiSquareStore;
	double kT_bins[kT_num_bins];
	
	for (int ii = 0; ii < kT_num_bins; ii++) {
		kT_bins[ii] = kT_min + ii * step;
	}
	
	for (int kTbin = 0; kTbin < kT_num_bins; kTbin++) {
		std::string filenames[6] = {
		"Correlation_Function_Dataset_Pionplus_" + std::to_string(kTbin) + ".data",
		"Correlation_Function_Dataset_Pionminus_" + std::to_string(kTbin) + ".data",
		"Correlation_Function_Dataset_Kaonplus_" + std::to_string(kTbin) + ".data",
		"Correlation_Function_Dataset_Kaonminus_" + std::to_string(kTbin) + ".data",
		"Correlation_Function_Dataset_Proton_" + std::to_string(kTbin) + ".data",
		"Correlation_Function_Dataset_Antiproton_" + std::to_string(kTbin) + ".data"
		};
		for (int fileIndex = 0; fileIndex < 6; fileIndex++) {
			TH3D *h3 = new TH3D("h3_bin", "Bin", q_num_bins, qmin, qmax, q_num_bins, qmin, qmax, q_num_bins, qmin, qmax);
			h3->Reset("ICE");
			std::ifstream file(filenames[fileIndex]);
			double qout_val, qside_val, qlong_val, Cq_val, err_val;
			while (file >> qout_val >> qside_val >> qlong_val >> Cq_val >> err_val) {
				h3->Fill(qout_val, qside_val, qlong_val, Cq_val);
				if (err_val > 0) {
    					h3->SetBinError(h3->FindBin(qout_val, qside_val, qlong_val), err_val);
				}
			}
			file.close();
			
			TF3 *fitFunc = new TF3("fitFunc", model_func_3D, qmin, qmax, qmin, qmax, qmin, qmax, 5);
			fitFunc->SetParLimits(1, Rout_min, Rout_max);
			fitFunc->SetParLimits(2, Rside_min, Rside_max);
			fitFunc->SetParLimits(3, Rlong_min, Rlong_max);
			fitFunc->FixParameter(0, 1.0);
			fitFunc->SetParameter(0, lambda_ini);
			fitFunc->SetParameter(1, Rout_ini);
			fitFunc->SetParameter(2, Rside_ini);
			fitFunc->SetParameter(3, Rlong_ini);
			
			h3->Fit(fitFunc, "QRB");
			
			double Rout = fitFunc->GetParameter(1);
			double Rside = fitFunc->GetParameter(2);
			double Rlong = fitFunc->GetParameter(3);
			
			int ndf = h3->GetNbinsX()*h3->GetNbinsY()*h3->GetNbinsZ() - fitFunc->GetNumberFreeParameters();
			double reducedChiSquare = fitFunc->GetChisquare() / ndf;
			
			std::cout << fitFunc->GetChisquare() << "      " << ndf << "      " << reducedChiSquare << std::endl ;

			// Write to file
			outputFile << kT_bins[kTbin] << "," << Rout << "," << Rside << "," << Rlong << "," << fileIndex << "\n";

			if (fileIndex == 0) {
			    Rout_storeval_1.push_back(Rout);
			    Rside_storeval_1.push_back(Rside);
			    Rlong_storeval_1.push_back(Rlong);
			} else if (fileIndex == 1) {
			    Rout_storeval_2.push_back(Rout);
			    Rside_storeval_2.push_back(Rside);
			    Rlong_storeval_2.push_back(Rlong);
			} else if (fileIndex == 2) {
			    Rout_storeval_3.push_back(Rout);
			    Rside_storeval_3.push_back(Rside);
			    Rlong_storeval_3.push_back(Rlong);
			} else if (fileIndex == 3) {
			    Rout_storeval_4.push_back(Rout);
			    Rside_storeval_4.push_back(Rside);
			    Rlong_storeval_4.push_back(Rlong);
			} else if (fileIndex == 4) {
			    Rout_storeval_5.push_back(Rout);
			    Rside_storeval_5.push_back(Rside);
			    Rlong_storeval_5.push_back(Rlong);
			} else {
			    Rout_storeval_6.push_back(Rout);
			    Rside_storeval_6.push_back(Rside);
			    Rlong_storeval_6.push_back(Rlong);
			}
			
			delete fitFunc;
			delete h3;
		}
	}
	outputFile.close();
	
	TCanvas *c1 = new TCanvas("c1", "Radii_vs_kT_plots", 1200, 400);
	c1->Divide(3, 1, 0.001, 0.001);
	TGraph *graphs[18] = {
	new TGraph(kT_num_bins, kT_bins, Rout_storeval_1.data()),
	new TGraph(kT_num_bins, kT_bins, Rside_storeval_1.data()),
	new TGraph(kT_num_bins, kT_bins, Rlong_storeval_1.data()),
	new TGraph(kT_num_bins, kT_bins, Rout_storeval_2.data()),
	new TGraph(kT_num_bins, kT_bins, Rside_storeval_2.data()),
	new TGraph(kT_num_bins, kT_bins, Rlong_storeval_2.data()),
	new TGraph(kT_num_bins, kT_bins, Rout_storeval_3.data()),
	new TGraph(kT_num_bins, kT_bins, Rside_storeval_3.data()),
	new TGraph(kT_num_bins, kT_bins, Rlong_storeval_3.data()),
	new TGraph(kT_num_bins, kT_bins, Rout_storeval_4.data()),
	new TGraph(kT_num_bins, kT_bins, Rside_storeval_4.data()),
	new TGraph(kT_num_bins, kT_bins, Rlong_storeval_4.data()),
	new TGraph(kT_num_bins, kT_bins, Rout_storeval_5.data()),
	new TGraph(kT_num_bins, kT_bins, Rside_storeval_5.data()),
	new TGraph(kT_num_bins, kT_bins, Rlong_storeval_5.data()),
	new TGraph(kT_num_bins, kT_bins, Rout_storeval_6.data()),
	new TGraph(kT_num_bins, kT_bins, Rside_storeval_6.data()),
	new TGraph(kT_num_bins, kT_bins, Rlong_storeval_6.data())
	};
	int colors[6] = {kRed, kViolet, kBlue, kBlue-6, kGreen, kGreen+4};
	const char *titles[3] = {"R_{out}", "R_{side}", "R_{long}"};
	
	graphs[0]->GetXaxis()->SetTitle("k_{T} [GeV/c]");
	graphs[1]->GetXaxis()->SetTitle("k_{T} [GeV/c]");
	graphs[2]->GetXaxis()->SetTitle("k_{T} [GeV/c]");
	
	graphs[0]->GetYaxis()->SetTitle("R_{out} [fm]");
	graphs[1]->GetYaxis()->SetTitle("R_{side} [fm]");
	graphs[2]->GetYaxis()->SetTitle("R_{long} [fm]");
	
	graphs[0]->GetXaxis()->SetLimits(0.01, 1.51);
	graphs[1]->GetXaxis()->SetLimits(0.01, 1.51);
	graphs[2]->GetXaxis()->SetLimits(0.01, 1.51);
	
	for (int i = 0; i < 3; i++) {
		c1->cd(i + 1);
		
		graphs[i]->SetTitle((std::string(titles[i]) + " vs k_{T}").c_str());
		
		graphs[i]->SetMarkerColor(colors[0]);
		graphs[i]->SetLineColor(colors[0]);
		graphs[i]->SetLineStyle(1);
		graphs[i]->SetLineWidth(1);

		graphs[i + 3]->SetMarkerColor(colors[1]);
		graphs[i + 3]->SetLineColor(colors[1]);
		graphs[i + 3]->SetLineStyle(2);
		graphs[i + 3]->SetLineWidth(1);
		
		graphs[i + 6]->SetMarkerColor(colors[2]);
		graphs[i + 6]->SetLineColor(colors[2]);
		graphs[i + 6]->SetLineStyle(1);
		graphs[i + 6]->SetLineWidth(1);
		
		graphs[i + 9]->SetMarkerColor(colors[3]);
		graphs[i + 9]->SetLineColor(colors[3]);
		graphs[i + 9]->SetLineStyle(2);
		graphs[i + 9]->SetLineWidth(1);
		
		graphs[i + 12]->SetMarkerColor(colors[4]);
		graphs[i + 12]->SetLineColor(colors[4]);
		graphs[i + 12]->SetLineStyle(1);
		graphs[i + 12]->SetLineWidth(1);
		
		graphs[i + 15]->SetMarkerColor(colors[5]);
		graphs[i + 15]->SetLineColor(colors[5]);
		graphs[i + 15]->SetLineStyle(2);
		graphs[i + 15]->SetLineWidth(1);

		if (i==0 || i== 1) {
			graphs[i]->SetMinimum(3.0);
			graphs[i]->SetMaximum(7.0);
			graphs[i+3]->SetMinimum(3.0);
			graphs[i+3]->SetMaximum(7.0);
			graphs[i+6]->SetMinimum(3.0);
			graphs[i+6]->SetMaximum(7.0);
			graphs[i+9]->SetMinimum(3.0);
			graphs[i+9]->SetMaximum(7.0);
			graphs[i+12]->SetMinimum(3.0);
			graphs[i+12]->SetMaximum(7.0);
			graphs[i+15]->SetMinimum(3.0);
			graphs[i+15]->SetMaximum(7.0);
		} else {
			graphs[i]->SetMinimum(3.0);
			graphs[i]->SetMaximum(11.0);
			graphs[i+3]->SetMinimum(3.0);
			graphs[i+3]->SetMaximum(11.0);
			graphs[i+6]->SetMinimum(3.0);
			graphs[i+6]->SetMaximum(11.0);
			graphs[i+9]->SetMinimum(3.0);
			graphs[i+9]->SetMaximum(11.0);
			graphs[i+12]->SetMinimum(3.0);
			graphs[i+12]->SetMaximum(11.0);
			graphs[i+15]->SetMinimum(3.0);
			graphs[i+15]->SetMaximum(11.0);
		}
			
		graphs[i]->Draw("APC");
		graphs[i+3]->Draw("PC same"); 
		graphs[i+6]->Draw("PC same"); 
		graphs[i+9]->Draw("PC same"); 
		graphs[i+12]->Draw("PC same");
		graphs[i+15]->Draw("PC same"); 
		
		TLegend *leg = new TLegend(0.7, 0.5, 0.85, 0.85);
		leg->AddEntry(graphs[i], "  #pi^{+}-#pi^{+}", "l");
		leg->AddEntry(graphs[i+3], "  #pi^{-}-#pi^{-}", "l");
		leg->AddEntry(graphs[i+6], "  #kappa^{+}-#kappa^{+}", "l");
		leg->AddEntry(graphs[i+9], "  #kappa^{-}-#kappa^{-}", "l");
		leg->AddEntry(graphs[i+12], "  p-p", "l");
		leg->AddEntry(graphs[i+15], "  #bar{p}-#bar{p}", "l");
		leg->SetBorderSize(0);
		leg->SetTextSize(0.06);
		leg->Draw();
	}
	c1->SaveAs("Radii_vs_kT_plots_All_.eps");
	delete c1;
	for (TGraph *g : graphs) delete g;
	return 0;
}
