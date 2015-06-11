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
#include "RooAddPdf.h"
#include "RooArgList.h"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
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
 
 /*
 //gives me the mass histogram under a cut applied 
TH1F* Mass(TString variable, Double_t xMin , Double_t xMax , const Char_t* Cut_or_Nah)
{
	cout << "Getting mass" << endl;
	// Start by looping over the Photon + Jet 
	TFile *inclusive = new TFile("/home/bbachu/cms/cmssw/040/CMSSW_7_4_0/src/Tag_and_Probe/skimmed.root");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny("events");
	Int_t pair;

	Float_t tag_px[64] ;
	Float_t tag_py[64] ;
	Float_t tag_pz[64] ;
	Float_t tag_energy[64] ;

	Float_t probe_px[64] ;
	Float_t probe_py[64] ;
	Float_t probe_pz[64] ;
	Float_t probe_energy[64] ;
	UInt_t nPairs;

	tree->SetBranchAddress("tag.px" , &tag_px);
	tree->SetBranchAddress("tag.py" , &tag_py);
	tree->SetBranchAddress("tag.pz", &tag_pz);
	tree->SetBranchAddress("tag.energy", &tag_energy);

	tree->SetBranchAddress("probe.px", &probe_px );
	tree->SetBranchAddress("probe.py", &probe_py);
	tree->SetBranchAddress("probe.pz", &probe_pz);
	tree->SetBranchAddress("probe.energy", &probe_energy);
	tree->SetBranchAddress("nPairs",&nPairs);

	//set up the variable to cut over
	Float_t variable;
	tree->SetBranchAddress(variable , &variable);
 
	//create hist to store the mass
	TH1F *h_mass = new TH1F( "h_mass","Invariant Mass", 100, 0, 200);
	//read all entries and fill the hist
	Long64_t nEvents = (Int_t) tree->GetEntries();
	Double_t mass;

	for (Long64_t entry = 0; entry < nEvents; entry++ )
	{
		tree->GetEntry(entry);
		//lopp through pairs
		for (Int_t ipair = 0 ; ipair < nPairs ; ipair++)
		{
			//check if CUT was specified, if it was start applying cut over range
			if ( Cut_or_Nah == "Cut") and !(  (xMin < variable) && (variable < xMax) ) continue;
			//construct the 4 vector for the tag
			TLorentzVector *v_tag = new TLorentzVector();
			v_tag->SetPxPyPzE(tag_px[ipair], tag_py[ipair], tag_pz[ipair], tag_energy[ipair]);
			//construct the 4 vector for the probe
			TLorentzVector *v_probe = new TLorentzVector();
			v_probe->SetPxPyPzE(probe_px[ipair], probe_py[ipair], probe_pz[ipair], probe_energy[ipair]);
			//construct mass
			mass = 	(*v_probe+ *v_tag ).M() ;
			h_mass->Fill(mass);
		}
	}
	cout << "exit loop" << endl;

	//TCanvas *c1 = new TCanvas();
	//h_mass->Draw();
	cout << "Got all the mass" << endl;
return h_mass;
}

*/


//gives me the mass histogram under a cut applied 
template<typename T> TH1F* Mass(TString variable, Double_t xMin , Double_t xMax , const Char_t *Cut_or_Nah , Int_t ElectronVeto_or_nah)
{
	cout << "Getting mass" << endl;
	// Start by looping over the Photon + Jet 
	TFile *inclusive = new TFile("/home/bbachu/cms/cmssw/040/CMSSW_7_4_0/src/Tag_and_Probe/skimmed.root");
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
	T x;
	tree->SetBranchAddress(variable, &x);
	
	//create hist to store the mass
	TH1F *h_mass = new TH1F( "h_mass","Invariant Mass", 100, 0, 200);
	//read all entries and fill the hist
	Long64_t nEvents = (Int_t) tree->GetEntries();

	for (Long64_t entry = 0; entry < nEvents; entry++ )
	{
		tree->GetEntry(entry);		
		//check if we should apply elecron veto
		if( ElectronVeto_or_nah != electronVetoApplied ) continue;
			//check if CUT was specified, if it was start applying cut over range
			if ( Cut_or_Nah == "Cut")
			{
				if( !((xMin < x) && (x < xMax)) ) continue;
				h_mass->Fill(mass);
			} 
			else
			{
				h_mass->Fill(mass);
			}
		
	}
	cout << "exit loop" << endl;

	//TCanvas *cMass = new TCanvas();
	//h_mass->Draw();
	cout << "Got all the mass" << endl;
return h_mass;
}

//fit for Z mass
void Fit_for_Z(TH1F* h)
{
	//TCanvas *c2 = new TCanvas();
	//observable for the pdf
	RooRealVar *Z_mass = new RooRealVar("Z_mass" , "Z_mass" , 0, 300);
	
	//create the Breit Wigner function SIGNAL
	RooRealVar *m = new RooRealVar("m", "m", 91.1876 , 87 , 95);
	//m->setConstant(kTrue);
	RooRealVar *g = new RooRealVar("g", "g", 200, 0, 6000);
	RooBreitWigner *Z_Peak = new RooBreitWigner("Z_Peak", "Z_Peak", *Z_mass, *m , *g);

	//create the Exp for the BACKGROUND
	RooRealVar *a = new RooRealVar("a","a",0 ,0 ,-3 );
	RooExponential *bkg = new RooExponential("bkg", "bkg", *Z_mass, *a );

	RooRealVar *nsig = new RooRealVar("nsig", "signal fraction", 25000 , 0.  , 26495.);
	RooRealVar *nbkg = new RooRealVar("nbkg", "background fraction" , 5000 ,  0.  , 26495.);

	RooAddPdf *model = (RooAddPdf*) new RooAddPdf("model" , "model" , RooArgList(*Z_Peak, *bkg) , RooArgList(*nsig, *nbkg));

	//make the histogram for the fitting
	RooDataHist *data_Z_mass = (RooDataHist*) new RooDataHist("data", "dataset with mass", *Z_mass, h);

	//fit the model to the data
	model->fitTo(*data_Z_mass , Extended(kTRUE));

	//plot pdf and data overlaid
	RooPlot *xframe = Z_mass->frame();
	data_Z_mass->plotOn(xframe);
	model->plotOn(xframe);
	model->plotOn(xframe, Components(*Z_Peak), LineStyle(kDashed), LineColor(kRed));
	model->plotOn(xframe, Components(*bkg), LineStyle(kDotted), LineColor(kGreen));
	//xframe->Draw();

	//what is the structure of my composite model
	model->printCompactTree();

	RooArgSet *params = model->getVariables();
	//get the number of entries in the BreitWeigner Distribution
	RooRealVar *integral_BW = (RooRealVar*) params->find("nsig");
	Double_t Integral_BW = integral_BW->getValV();
	cout <<"Integral over BreitWeigner = " << Integral_BW << endl;
}

//this model uses the extended likelihood formalism and a composite pdf
//the composite pdf is composed of a Decaying exponential and BreitWeigner Distribution
RooAbsPdf* build_model(RooRealVar& mass)
{
	//create the Breit Wigner function SIGNAL
	RooRealVar *m = new RooRealVar("m", "m", 91 , 87 , 95);
	//m->setConstant(kTrue);
	RooRealVar *g = new RooRealVar("g", "g", 200, 0, 6000);
	RooBreitWigner *Z_Peak = new RooBreitWigner("Z_Peak", "Z_Peak", mass, *m , *g);

	//create the Exp for the BACKGROUND
	RooRealVar *a = new RooRealVar("a","a",0 ,0 ,-3 );
	RooExponential *bkg = new RooExponential("bkg", "bkg", mass, *a );

	RooRealVar *nsig = new RooRealVar("nsig", "signal fraction", 25000 , 0.  , 26495.);
	RooRealVar *nbkg = new RooRealVar("nbkg", "background fraction" , 5000 ,  0.  , 26495.);

	RooAddPdf *model = (RooAddPdf*) new RooAddPdf("model" , "model" , RooArgList(*Z_Peak, *bkg) , RooArgList(*nsig, *nbkg));

return model;

}

void Print_chiSquare(RooPlot *frame, const char* histname , const char* pdfname)
{	
	Double_t ChiSquare;
	ChiSquare =  frame->chiSquare( pdfname , histname , 8);
	cout << "ChiSquare result from " << pdfname << " = "  << ChiSquare << endl;
}

//this model will use a convolution of the BreitWeigner Distribution with the Crystal Ball Distribution as the fit for the Z peak
RooAddPdf* build_model_2(RooRealVar& mass)
{
	cout << "Building Model to work with Detector Resolution" << endl ;
	//PHYSICS MODEL
	//create the Breit Wigner function for the Z peak
	RooRealVar *m = new RooRealVar("m", "m", 91 , 87 , 95);
	RooRealVar *g = new RooRealVar("g", "g", 200, 0, 6000);
	RooBreitWigner *Z_Peak_model = new RooBreitWigner("Z_Peak_model", "Z Peak Model", mass, *m , *g);

	//DETECTOR RESOLUTION MODEL
	//create Crystal Ball function for the the Z peak
	RooRealVar *m2 = new RooRealVar("m2", "m2", 0 , -0.00009, 0.00009);
	//m2->setConstant(kTRUE);
	RooRealVar *s = new RooRealVar("s", "s",1 , -100, 100 );
	RooRealVar *a2 = new RooRealVar("a2", "a2" , 1000 , 0 , 10000 );
	//a2->setConstant(kTRUE);
	RooRealVar *n = new RooRealVar("n", "n", 1, 0, 30000 );
	RooCBShape *Detec_Res = new RooCBShape("Detec_Res", "Detector Resolution", mass, *m2, *s, *a2, *n);

	cout << "Creating convolution of the physics MODEL and Resolution" << endl;
	//construction the convolution
	RooFFTConvPdf *MODELxRES = new RooFFTConvPdf("MODELxRES" , " Model: BreitWeigner (x) Resolution: Crystal Ball" , mass , *Detec_Res ,*Z_Peak_model  );

	//create the Exp for the BACKGROUND
	RooRealVar *a = new RooRealVar("a","a",0 ,0 ,-3 );
	RooExponential *bkg = new RooExponential("bkg", "bkg", mass, *a );

	RooRealVar *nsig = new RooRealVar("nsig", "signal fraction", 25000, 0.  , 26495);
	RooRealVar *nbkg = new RooRealVar("nbkg", "background fraction" , 5000 ,  0.  , 26495.);

	RooAddPdf *model = (RooAddPdf*) new RooAddPdf("model" , "model" , RooArgList(*MODELxRES, *bkg) , RooArgList(*nsig , *nbkg ));

return model;
}

Double_t Fit_model_to_data_and_integrate(TH1F* h)
{
	
	//observable for the pdf
	//RooRealVar *Z_mass = new RooRealVar("Z_mass" , "Z_mass" , 0, 300);
	RooRealVar mass("mass" , "mass" , 0, 300);
	RooAbsPdf *model = build_model_2(mass);
	
	//make the histogram for the fitting
	RooDataHist *data_mass = (RooDataHist*) new RooDataHist("data", "dataset with mass", mass, h);
	cout << "Fitting..." << endl;
	//fit the model to the data
	model->fitTo(*data_mass , Extended(kTRUE) , Minimizer("Minuit2", "Migrad"));

	//TCanvas *c3 = new TCanvas();
	//plot pdf and data overlaid
	RooPlot *xframe = mass.frame();
	data_mass->plotOn(xframe);
	model->plotOn(xframe, Name("model"));
	//model->plotOn(xframe, Components("*Z_Peak_model"), LineStyle(kDashed), LineColor(kRed));
	//model->plotOn(xframe, Components("*Detec_Res"), LineStyle(kDashed), LineColor(kYellow));
	//model->plotOn(xframe, Components("*bkg"), LineStyle(kDotted), LineColor(kGreen));
	//xframe->Draw();

	//what is the structure of my composite model
	model->printCompactTree();

	RooArgSet *params = model->getVariables();
	//get the number of entries in the BreitWeigner Distribution
	RooRealVar *integral_BWxCB = (RooRealVar*) params->find("nsig");
	Double_t Integral_BWxCB = integral_BWxCB->getValV();
	cout <<"Integral over BreitWeigner = " << Integral_BWxCB << endl;

	Print_chiSquare( xframe , "h_data" , "model");

	xframe->Print("v");

return Integral_BWxCB;
}

//get the ratio of the integrals when the electron veto is applied versus when it is not
Double_t Get_Fake_Rate(TString variable, Double_t x1 ,Double_t x2 )
{
	TH1F* h_cut_on_probe;
	TH1F* h_cut_on_probe_electron_ID;
	if (variable=="nVertices"){
		//produce histogram with cut on probe
		h_cut_on_probe = (TH1F*) Mass<UInt_t>( variable,  x1 , x2 , "Cut" , 0 );
		//produce histogram with cut on probe and also electron ID
		h_cut_on_probe_electron_ID = (TH1F*) Mass<UInt_t>( variable , x1, x2, "Cut", 1 );
	}
	else {
		//produce histogram with cut on probe
		h_cut_on_probe = (TH1F*) Mass<Float_t>( variable,  x1 , x2 , "Cut" , 0 );
		//produce histogram with cut on probe and also electron ID
		h_cut_on_probe_electron_ID = (TH1F*) Mass<Float_t>( variable , x1, x2, "Cut", 1 );
	}
	//now I have two histograms, one with Electron ID applied and one without

	//get the integrals required from both histograms
	Double_t integral_phoID = Fit_model_to_data_and_integrate(h_cut_on_probe);
	Double_t integral_phoID_and_electronID = Fit_model_to_data_and_integrate(h_cut_on_probe_electron_ID);

	cout << "Calculating fake rate for " << variable << endl;
	//get fake rate
	Double_t fake_rate = (integral_phoID_and_electronID / integral_phoID) ;
return fake_rate;
}

//Plot fake rate as a function of variable
TH1F* fake_rate_phase_space(TString variable , Int_t nbin , Double_t xMin , Double_t xMax)
{
	//TCanvas *c = new TCanvas();
	//calculate size of increments
	Double_t increment = (xMax - xMin) / nbin ;
	//create a histogram to represent the fake rate versus variable
	TH1F* hEfficiency = (TH1F*) new TH1F(variable, "Efficiency as a function of " + variable, nbin , xMin , xMax);
	//loop over bins of the variable in question to plot efficiency 
	for( Int_t i = 0 ; i < nbin-1 ; i++)
	{
		Double_t fake_rate = Get_Fake_Rate( variable, xMin + (i*increment) , increment + (i*increment) );
		hEfficiency->SetBinContent( i+1 , fake_rate);
	}
	//hEfficiency->Draw();
	cout << "Plotting Efficiency versus " << variable << endl;

return hEfficiency;
}

//make a good draw function to show the fake rate versus variable
void Draw_phase_space( TString variable , Int_t nbins , Double_t xMin , Double_t xMax )
{
	cout << "Beginning phase space calculation for " << variable << endl;
	TCanvas *c = new TCanvas();
	//create the histogram of the Efficiency versus variable
	TH1F* h_phase_space = (TH1F*) fake_rate_phase_space( variable , nbins , xMin , xMax);
	h_phase_space->GetXaxis()->SetTitle(variable.Data());
	h_phase_space->GetYaxis()->SetTitle("Fake Rate");
	//const Char_t* title =  ;
	h_phase_space->SetTitle( "Fake Rate as a function of " + variable);
	h_phase_space->SetStats(0);
	cout << "Finished phase space for " << variable << endl;
	h_phase_space->Draw();
}


//main main function 
void e_gamma_tandp()
{
	cout <<"Starting main function" << endl;
	//TH1F *h_mass = (TH1F*) Mass();
	//Fit_for_Z(h_mass);
	//Fit_model_to_data(h_mass);

	//reproduce the first histogam and fit we made
	// TCanvas *initial = new TCanvas();
	// cout << "Reproducing first fit as a check" << endl;
	// TH1F *h_mass = (TH1F*) Mass(" ", 0 , 0 , "");
	// h_mass->Draw();

	//make phase space histograms
	cout<< "About to draw phase space diagrams" << endl;
	Draw_phase_space( "p_T" , 10 , 0 , 250 );
	Draw_phase_space("eta", 10, -2.5, 2.5 );
	Draw_phase_space("nVertices" , 10 , 0 , 30);
	cout <<"End main function" << endl;

}