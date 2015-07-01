/*gSystem->AddIncludePath(TString("-l/cvmfs/cms.cern.ch/")+TString(gSystem->Getenv("SCRAM_ARCH"))+
			       TString("/lcg/roofit/5.32.00/include/"));*/
//#include "RooRealVar.h"
//#include "RooDataSet.h"
//#include "RooGaussian.h"
//#include "RooExponential.h"
#include "TCanvas.h"
// #include "RooPlot.h"
#include "TAxis.h"
#include "TCut.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TImage.h"
#include "TROOT.h"
#include "fstream"
#include "string"
#include "sstream"
#include "iostream"
#include "iomanip"
#include "TCut.h"
#include "vector"
#include "TLegend.h"
#include "TMath.h"
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
//#include "RooAddPdf.h"
//#include "vector.h"
//#include "My_Functions.h"
//#include "Draw_Histograms.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "RooGlobalFunc.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooArgList.h"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace RooFit ;
TString outPath = "~/local/Dropbox/mit/figs/egamma_v5";
TFile * fout = new TFile(outPath+"/output.root","RECREATE");

//gives me the mass histogram under a cut applied 
TH1F* Mass(TString variable, Double_t xMin , Double_t xMax , const Char_t *Cut_or_Nah , Int_t ElectronVeto_or_nah)
{
	cout << "Getting mass" << endl;
	// Start by looping over the Photon + Jet 
	TFile *inclusive = new TFile("/home/snarayan/cms/hist/egamma_v5/t2mit/filefi/032/flattened.root");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny("events");

	Float_t mass ;
	// Float_t eta ;
	// Float_t p_T ;
	Int_t electronVetoApplied ;

	tree->SetBranchAddress("mass" , &mass);
	tree->SetBranchAddress("electronVetoApplied" , &electronVetoApplied);
	// tree->SetBranchAddress("eta" , &eta);
	// tree->SetBranchAddress("p_T", &p_T);
	// tree->SetBranchAddress("electronVetoApplied", &electronVetoApplied);

	//set up the variable to cut over
	// T x;
	// tree->SetBranchAddress(variable, &x);
	
	//create hist to store the mass
	// TH1F *h_mass = new TH1F( "h_mass","Invariant Mass", 100, 0, 200);
	TH1F *h_mass = new TH1F("h_mass","h_mass",50,60,120);
	//read all entries and fill the hist
	if (ElectronVeto_or_nah==0){
		if (Cut_or_Nah == "Cut")
			tree->Draw("mass>>h_mass",TString::Format("%f < %s && %s < %f && electronVetoApplied==0 && mass > 0.01",xMin,variable.Data(),variable.Data(),xMax),"goff");
		else
			tree->Draw("mass>>h_mass","electronVetoApplied==0 && mass > 0.01","goff");
	} else {
		if (Cut_or_Nah == "Cut")
			tree->Draw("mass>>h_mass",TString::Format("%f < %s && %s < %f && electronVetoApplied==1 && mass > 0.01",xMin,variable.Data(),variable.Data(),xMax),"goff");
		else
			tree->Draw("mass>>h_mass","electronVetoApplied==1 && mass > 0.01","goff");
	}
	
	cout << "Got all the mass" << endl;
return h_mass;
}

//this model will use a convolution of the BreitWeigner Distribution with the Crystal Ball Distribution as the fit for the Z peak
RooAddPdf* build_model_2(RooRealVar& mass , Double_t nEntries)
{
	cout << "Building Model to work with Detector Resolution" << endl ;
	//PHYSICS MODEL
	//create the Breit Wigner function for the Z peak
	RooRealVar *m = new RooRealVar("m", "m", 91. , 88. , 94.);
	RooRealVar *g = new RooRealVar("g", "g", 10., 0., 100.);
	RooBreitWigner *Z_Peak_model = new RooBreitWigner("Z_Peak_model", "Z Peak Model", mass, *m , *g);

	//DETECTOR RESOLUTION MODEL
	//create Crystal Ball function for the the Z peak
	// RooRealVar *m2 = new RooRealVar("m2", "m2", 0 , -0.00009, 0.00009);
	// //m2->setConstant(kTRUE);
	// RooRealVar *s = new RooRealVar("s", "s",1 , -100, 100 );
	// RooRealVar *a2 = new RooRealVar("a2", "a2" , 1000 , 0 , 10000 );
	// //a2->setConstant(kTRUE);
	// RooRealVar *n = new RooRealVar("n", "n", 1, 0, 30000 );
	// RooCBShape *Detec_Res = new RooCBShape("Detec_Res", "Detector Resolution", mass, *m2, *s, *a2, *n);

	//DETECTOR RESOLUTION MODEL with gaussian (CB is confusing)
	RooRealVar *mean = new RooRealVar("mean" , "mean" , 0. , -2 , 2);
	RooRealVar *sigma = new RooRealVar("sigma", "sigma" , 5. , 0., 30.);
	RooGaussian *Gaussian_resolution = new RooGaussian("Gaussian_resolution", "Detector Resolution", mass , *mean , *sigma);

	cout << "Creating convolution of the physics MODEL and Resolution" << endl;
	//construction the convolution
	RooFFTConvPdf *MODELxRES = new RooFFTConvPdf("MODELxRES" , " Model: BreitWeigner (x) Resolution: Gaussian" , mass , *Gaussian_resolution ,*Z_Peak_model  );

	//create the Exp for the BACKGROUND
	RooRealVar *a = new RooRealVar("a","a",-0.001 ,-3. ,0. );
	RooExponential *bkg = new RooExponential("bkg", "bkg", mass, *a );

	RooRealVar *a0 = new RooRealVar("a0","a0",1,-100,100);
	RooRealVar *a1 = new RooRealVar("a1","a1",1,-100,100);
	RooRealVar *a2 = new RooRealVar("a2","a2",1,-100,100);
	RooRealVar *a3 = new RooRealVar("a3","a3",1,-100,100);
	RooRealVar *a4 = new RooRealVar("a4","a4",1,-100,100);
	RooPolynomial *poly = new RooPolynomial("poly","poly",mass,RooArgList(*a0,*a1,*a2, *a3),0);

	RooRealVar *nsig = new RooRealVar("nsig", "signal fraction", nEntries*0.8, nEntries*0.1  , nEntries);
	RooRealVar *nbkg = new RooRealVar("nbkg", "background fraction" , nEntries*0.2 ,  nEntries*.001  , nEntries);

	RooRealVar *sigFrac = new RooRealVar("sigFrac","signal fraction",.5,0,1);

	RooAddPdf *model = (RooAddPdf*) new RooAddPdf("model" , "model" , RooArgList(*MODELxRES, *bkg) , RooArgList(*nsig , *nbkg  ));
	// RooAddPdf *model = (RooAddPdf*) new RooAddPdf("model" , "model" , RooArgList(*Z_Peak_model, *bkg) , RooArgList(*nsig , *nbkg  ));

return model;
}

RooRealVar * Fit_model_to_data_and_integrate(TH1F* h , TString variable , Double_t x1 , Double_t x2)
{
	
	//observable for the pdf
	//RooRealVar *Z_mass = new RooRealVar("Z_mass" , "Z_mass" , 0, 300);
	RooRealVar mass("mass" , "mass" , 0, 300);
	Double_t nEntries = h->Integral();
	RooAbsPdf *model = build_model_2(mass , nEntries);
	
	//make the histogram for the fitting
	RooDataHist *data_mass = (RooDataHist*) new RooDataHist("data", "dataset with mass", mass, h);
	cout << "Fitting..." << endl;
	//fit the model to the data
	model->fitTo(*data_mass , Extended(kTRUE) , Minimizer("Minuit2", "Migrad"));
	printf("nEntries: %.0f\n",nEntries);
	// model->fitTo(*data_mass , Extended(kTRUE) , Minimizer("Minuit2", "Migrad"));

	// std::stringstream s;
	// s << "Mass Fit under cut" << x1 << " < " << variable << " < " << x2 ;
	// const char* str = s.str().c_str();

	TString s ="" ;
	s+= x1 ;
	s+= " < " ;
	s+= variable ;
	s+= " < " ;
	s+= x2 ;
	TCanvas *c3 = new TCanvas();
	//plot pdf and data overlaid
	RooPlot *xframe = mass.frame(Title(s));
	//xframe->SetNameTitle(  "LALALA" , "LALALA" );
	data_mass->plotOn(xframe);
	model->plotOn(xframe, Name("model"));
	//model->plotOn(xframe, Components("*Z_Peak_model"), LineStyle(kDashed), LineColor(kRed));
	//model->plotOn(xframe, Components("*Detec_Res"), LineStyle(kDashed), LineColor(kYellow));
	//model->plotOn(xframe, Components("*bkg"), LineStyle(kDotted), LineColor(kGreen));
	xframe->Draw("goff");
	// fout->cd();
	c3->SaveAs(outPath+TString::Format("/%s_%.3f_%.3f.png",variable.Data(),x1,x2));

	//what is the structure of my composite model
	model->printCompactTree();

	RooArgSet *params = model->getVariables();
	//get the number of entries in the BreitWeigner Distribution
	RooRealVar *integral_BWxCB = (RooRealVar*) params->find("nsig");
	Double_t Integral_BWxCB = integral_BWxCB->getValV();
	cout <<"Integral over BreitWeigner = " << Integral_BWxCB << endl;

	//Print_chiSquare( xframe , "h_data" , "model");
	xframe->Print("v");

	// // incase we decide to do integrals instead
	// RooAbsReal* integral_over_CBxBW = model->createIntegral(*mass );
	// Double_t Integral_BWxCB = integral_over_CBxBW->getVal();


return integral_BWxCB;
}

//get the ratio of the integrals when the electron veto is applied versus when it is not
std::vector<Double_t> Get_Fake_Rate(TString variable, Double_t x1 ,Double_t x2 )
{
	TH1F* h_cut_on_probe;
	TH1F* h_cut_on_probe_electron_ID;
	if (variable=="nVertices"){
		//produce histogram with cut on probe
		h_cut_on_probe = (TH1F*) Mass( variable,  x1 , x2 , "Cut" , 0 );
		//produce histogram with cut on probe and also electron ID
		h_cut_on_probe_electron_ID = (TH1F*) Mass( variable , x1, x2, "Cut", 1 );
	}
	else {	RooRealVar *sigFrac = new RooRealVar("sigFrac","signal fraction",.5,0,1);
		//produce histogram with cut on probe
		h_cut_on_probe = (TH1F*) Mass( variable,  x1 , x2 , "Cut" , 0 );
		//produce histogram with cut on probe and also electron ID
		h_cut_on_probe_electron_ID = (TH1F*) Mass( variable , x1, x2, "Cut", 1 );
	}
	//now I have two histograms, one with Electron ID applied and one without

	//get the integrals required from both histograms
	printf("===============\nSTARTING NO VETO\n==================\n");
	RooRealVar *integral_phoID = Fit_model_to_data_and_integrate(h_cut_on_probe ,  variable ,  x1 ,  x2);
	printf("===============\nSTARTING YES VETO\n==================\n");
	RooRealVar *integral_phoID_and_electronID = Fit_model_to_data_and_integrate(h_cut_on_probe_electron_ID ,  variable+TString("Veto") ,  x1 ,  x2);

	cout << "Calculating fake rate for " << variable << endl;
	//get fake rate
	Double_t N_sig_with_veto = integral_phoID_and_electronID->getValV() ;
	Double_t N_sig = integral_phoID->getValV();
	//compute fake rate
	Double_t fake_rate = ( N_sig_with_veto / N_sig ) ;
	Double_t N_sig_with_veto_Error = integral_phoID_and_electronID->getError();
	Double_t N_sig_Error = integral_phoID->getError();
	//compute error on fake rate
	Double_t Error_on_fake_rate = fake_rate * pow( pow( N_sig_Error/ N_sig , 2) + pow( N_sig_with_veto_Error/ N_sig_with_veto , 2 ) , 0.5);

	std::vector<Double_t> fake_rate_with_error;
	
	fake_rate_with_error.push_back(fake_rate);
	fake_rate_with_error.push_back(Error_on_fake_rate);

return fake_rate_with_error;
}

//Plot fake rate as a function of variable
std::vector<TH1F*> fake_rate_phase_space(TString variable , Int_t nbin , Double_t xMin , Double_t xMax)
{
	//TCanvas *c = new TCanvas();
	//calculate size of increments
	Double_t increment = (xMax - xMin) / nbin ;
	//create a histogram to represent the fake rate versus variable
	TH1F* hEfficiency = (TH1F*) new TH1F(variable, "Efficiency as a function of " + variable, nbin , xMin , xMax);
	TH1F* h_high_error = new TH1F(variable, "High error as a function of " + variable, nbin , xMin , xMax);
	TH1F* h_low_error = new TH1F(variable, "Low error as a function of " + variable, nbin , xMin , xMax);
	//loop over bins of the variable in question to plot efficiency 
	for( Int_t i = 0 ; i < nbin ; i++)
	{
		std::vector<Double_t> fake_rate_with_error = Get_Fake_Rate( variable, xMin + (i*increment) , xMin + increment + (i*increment) );
		Double_t fake_rate = fake_rate_with_error.at(0);
		Double_t fake_rate_plus_error = fake_rate_with_error.at(0) + fake_rate_with_error.at(1);
		Double_t fake_rate_minus_error = fake_rate_with_error.at(0) - fake_rate_with_error.at(1);
		// hEfficiency->SetBinContent( i+1 , fake_rate);
		hEfficiency->Fill(xMin+increment*(i+.5),fake_rate); cout << "fake_rate" <<fake_rate << endl;
		h_high_error->Fill(xMin+increment*(i+.5),fake_rate_plus_error); cout <<"fake_rate_plus_error"<< fake_rate_plus_error <<endl;
		h_low_error->Fill(xMin+increment*(i+.5),fake_rate_minus_error); cout<< "fake_rate_minus_error"<< fake_rate_minus_error<< endl;
	}
	//hEfficiency->Draw();
	cout << "Plotting Efficiency versus " << variable << endl;
	std::vector<TH1F*> v_of_hists;
	v_of_hists.push_back(hEfficiency);
	v_of_hists.push_back(h_high_error);
	v_of_hists.push_back(h_low_error);
return v_of_hists;
}

//make a good draw function to show the fake rate versus variable
void Draw_phase_space( TString variable , Int_t nbins , Double_t xMin , Double_t xMax )
{
	cout << "Beginning phase space calculation for " << variable << endl;
	TCanvas *c1 = new TCanvas();
	c1->cd();
	//create the histogram of the Efficiency versus variable
	std::vector<TH1F*> h_phase_space = fake_rate_phase_space( variable , nbins , xMin , xMax);
	h_phase_space.at(1)->SetLineColor(4);
	h_phase_space.at(2)->SetLineColor(2);
	h_phase_space.at(2)->Draw("hist");
	h_phase_space.at(1)->Draw("same hist");
	h_phase_space.at(0)->GetXaxis()->SetTitle(variable.Data());
	h_phase_space.at(0)->GetYaxis()->SetTitle("Fake Rate");
	//const Char_t* title =  ;
	h_phase_space.at(0)->SetLineColor(1);
	h_phase_space.at(0)->SetTitle( "Fake Rate as a function of " + variable);
	h_phase_space.at(0)->SetStats(0);
	cout << "Finished phase space for " << variable << endl;
	//h_phase_space.at(0)->Draw("same");
	fout->cd();
	h_phase_space.at(0)->Write(TString::Format("h_%s",variable.Data()));
	h_phase_space.at(1)->Write(TString::Format("h_%s_up",variable.Data()));
	h_phase_space.at(2)->Write(TString::Format("h_%s_down",variable.Data()));
}

//Plot fake rate as a function of variable
std::vector<TH1F*> fake_rate_phase_space(TString variable , Int_t nbin ,Double_t * bins)
{
	//TCanvas *c = new TCanvas();
	//calculate size of increments
	//create a histogram to represent the fake rate versus variable
	TH1F* hEfficiency =  new TH1F(variable, "Efficiency as a function of " + variable, nbin , bins);
	TH1F* h_high_error = new TH1F(variable, "High error as a function of " + variable, nbin , bins);
	TH1F* h_low_error = new TH1F(variable, "Low error as a function of " + variable, nbin , bins);
	//loop over bins of the variable in question to plot efficiency 
	for( Int_t i = 0 ; i < nbin ; i++)
	{
		std::vector<Double_t> fake_rate_with_error = Get_Fake_Rate( variable, bins[i] , bins[i+1] );
		Double_t fake_rate = fake_rate_with_error.at(0);
		Double_t fake_rate_plus_error = fake_rate_with_error.at(0) + fake_rate_with_error.at(1);
		Double_t fake_rate_minus_error = fake_rate_with_error.at(0) - fake_rate_with_error.at(1);
		// hEfficiency->SetBinContent( i+1 , fake_rate);
		Double_t fillVal = (bins[i]+bins[i+1])/2.;
		hEfficiency->Fill(fillVal,fake_rate); cout << "fake_rate" <<fake_rate << endl;
		h_high_error->Fill(fillVal,fake_rate_plus_error); cout <<"fake_rate_plus_error"<< fake_rate_plus_error <<endl;
		h_low_error->Fill(fillVal,fake_rate_minus_error); cout<< "fake_rate_minus_error"<< fake_rate_minus_error<< endl;
	}
	//hEfficiency->Draw();
	cout << "Plotting Efficiency versus " << variable << endl;
	std::vector<TH1F*> v_of_hists;
	v_of_hists.push_back(hEfficiency);
	v_of_hists.push_back(h_high_error);
	v_of_hists.push_back(h_low_error);
return v_of_hists;
}

//make a good draw function to show the fake rate versus variable
void Draw_phase_space( TString variable , Int_t nbins , Double_t *bins )
{
	cout << "Beginning phase space calculation for " << variable << endl;
	TCanvas *c1 = new TCanvas();
	c1->cd();
	//create the histogram of the Efficiency versus variable
	std::vector<TH1F*> h_phase_space = fake_rate_phase_space( variable , nbins , bins);
	h_phase_space.at(1)->SetLineColor(4);
	h_phase_space.at(2)->SetLineColor(2);
	h_phase_space.at(2)->Draw("hist");
	h_phase_space.at(1)->Draw("same hist");
	h_phase_space.at(0)->GetXaxis()->SetTitle(variable.Data());
	h_phase_space.at(0)->GetYaxis()->SetTitle("Fake Rate");
	//const Char_t* title =  ;
	h_phase_space.at(0)->SetLineColor(1);
	h_phase_space.at(0)->SetTitle( "Fake Rate as a function of " + variable);
	h_phase_space.at(0)->SetStats(0);
	cout << "Finished phase space for " << variable << endl;
	//h_phase_space.at(0)->Draw("same");
	fout->cd();
	h_phase_space.at(0)->Write(TString::Format("h_%s",variable.Data()));
	h_phase_space.at(1)->Write(TString::Format("h_%s_up",variable.Data()));
	h_phase_space.at(2)->Write(TString::Format("h_%s_down",variable.Data()));
}

//main main function 
void e_gamma_tandp()
{
	cout <<"Starting main function" << endl;

	//make phase space histograms
	cout<< "About to draw phase space diagrams" << endl;
	Int_t npTBins = 5;
	Double_t pTBins[] = {20,40,60,80,100,250};
	Draw_phase_space( "p_T" , npTBins, pTBins);
	Draw_phase_space("eta", 20, -2.5, 2.5 );
	// Draw_phase_space("eta" , 1 , 0 , .5);
	Draw_phase_space("nVertices" , 5 , 5,25);
	cout <<"End main function" << endl;

//	compute the total efficiency
	TH1F* h_total_without_electron_veto = (TH1F*) Mass( "p_T",  0 , 10 , "" , 0 );
	TH1F* h_total_with_electron_veto = (TH1F*) Mass( "p_T",  0 , 10 , " " , 1 );
	RooRealVar * nSigVeto = Fit_model_to_data_and_integrate(h_total_with_electron_veto ,"", 1 ,10 );
	RooRealVar * nSig = Fit_model_to_data_and_integrate(h_total_without_electron_veto , "" , 1 ,10);
	Double_t integralVeto = nSigVeto->getVal();
	Double_t integralVetoErr = nSigVeto->getError();
	Double_t integral = nSig->getVal();
	Double_t integralErr = nSig->getError();
	cout << "Integral with veto = " << integralVeto << " +/- " << integralVetoErr << endl;
	cout << "Integral without veto = " << integral << " +/- " << integralErr <<  endl;
	Double_t Total_Efficiency = integralVeto/integral ;	
	Double_t Total_EfficiencyErr = Total_Efficiency * TMath::Sqrt(pow(integralVetoErr/integralVeto,2) + pow(integralErr/integral,2));
	cout <<" Total_Efficiency Efficiency = " << Total_Efficiency << " +/- " << Total_EfficiencyErr << endl ;

	fout->Close();

}
