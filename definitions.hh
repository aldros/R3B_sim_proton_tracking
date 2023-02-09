#ifndef _DEFINITIONS_HH_
#define _DEFINITIONS_HH_

//Charge and mass of the simulatied ion
const double Qion = 8.0000000000;
const double Aion = 16.0000000000;

UInt_t Nstat = 50000;//events used for training

//------------- Set Multi Dimensional Fit Class for Tracking -----------
const UInt_t   nVars   	= 4;
const UInt_t   nReduct 	= 0;//numbers of vars to ignore in PCA 

const UInt_t   Power   	= 8;//
const Double_t Fraction	= 1.;//statistics fraction for training sample

const UInt_t   MDFMaxFunctions = 10000;
const UInt_t   MDFMaxStudy = 100000;
const UInt_t   MDFMaxTerms = 50;
const UInt_t   MDFPowerLimit = 1;
const Double_t MDFMinAngle = 0.1;
const Double_t MDFMaxAngle = 0.1;
const Double_t MDFMinError = 1e-9;

//Histograms for P/Q
//TH2F * h_track_vs_fit 	= new TH2F("h_track_vs_fit","h_track_vs_fit",1000,0,6,1000,0,6);
//TH1F * h_residual 	= new TH1F("h_residual","h_residual",1000,-0.1,0.1);
//TH1F * h_track		= new TH1F("h_track","h_track",1000,0,6);
//TH2F * h_residual_vs_track = new TH2F("h_residual_vs_track","h_residual_vs_track",
//        1000,0,6,1000,-0.1,0.1);

//histograms for Target X,Y
TH2F * h_track_vs_fit = new TH2F("h_track_vs_fit","h_track_vs_fit",1000,-10,10,1000,-10,10);
TH1F * h_residual 	= new TH1F("h_residual","h_residual",1000,-10,10);
TH1F * h_track 	= new TH1F("h_track","h_track",1000,-10,10);
TH2F * h_residual_vs_track	= new TH2F("h_residual_vs_track","h_residual_vs_track",1000,-10,10,1000,-10,10);

//histograms for Target TX,TY
//TH2F * h_track_vs_fit = new TH2F("h_track_vs_fit","h_track_vs_fit",1000,-0.1,0.1,1000,-0.1,0.1);
//TH1F * h_residual 	= new TH1F("h_residual","h_residual",1000,-0.05,0.05);
//TH1F * h_track	= new TH1F("h_track","h_track",1000,-0.1,0.1);
//TH2F * h_residual_vs_track	= new TH2F("h_residual_vs_track","h_residual_vs_track",1000,-0.2,0.2,1000,-0.05,0.05);

//TH2F * h_ZvsX	= new TH2F("h_ZvsX","h_ZvsX",2000,-200,800,1000,-300,100);
#endif
