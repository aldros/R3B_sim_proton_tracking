#include"libs.hh"
#include"definitions.hh"
#include"run_proton_tracking.hh"
#include <TMath.h>
#include "mass.h"
//#include"R3BMDFWrapper.h"
//#include"R3BMDFWrapper.cxx"
//class R3BMDFWrapper;

using namespace std;

Double_t fNonUniformity = 1.0; // Experimental non-uniformity parameter
Double_t fResolution = 0.01;    // Experimental resolution
Double_t fComponentRes = 0.005;  // Experimental resolution for Nf and Ns
Double_t fThreshold = 0.020;    // 0.000010 // Minimum energy requested to create a Cal

Double_t UNIT = 931.4940954;
Double_t Mp = 938.272;
Double_t Mfrag = 11.009305404*UNIT;
Double_t Mbeam = 12.*UNIT;

//Track variable is container for an event
struct Data
{
    Double_t edata[nVars];
    Double_t pdata[nVars - nReduct];
    Double_t value;
};

Double_t ElossToEnergy_PunchThrough(Double_t X){

  Double_t p0 = 4.188;
  Double_t p1 = -35.82;
  Double_t p2 = -0.6455;
  Double_t p3 = -1.73;
  
  Double_t res = TMath::Exp(p0+p1*X) + TMath::Exp(p2+p3*X);
  
  return res;
}

Double_t NUSmearing(Double_t inputEnergy){
    // Very simple preliminary scheme where the NU is introduced as a flat random
    // distribution with limits fNonUniformity (%) of the energy value.
    //
  gRandom = new TRandom3();
  gRandom->SetSeed(1);

  return gRandom->Uniform(inputEnergy - inputEnergy * fNonUniformity / 100,
			  inputEnergy + inputEnergy * fNonUniformity / 100);
}

Double_t ExpResSmearing(Double_t inputEnergy){
  // Smears the energy according to some Experimental Resolution distribution
  // Very simple  scheme where the Experimental Resolution
  // is introduced as a gaus random distribution with a width given by the
  // parameter fResolution(in % @ MeV). Scales according to 1/sqrt(E)
  gRandom = new TRandom3();
  gRandom->SetSeed(1);

  if (fResolution == 0){
    return inputEnergy;
  }else{
    Double_t randomIs = gRandom->Gaus(0, inputEnergy * fResolution / (sqrt(inputEnergy)));
    return inputEnergy + randomIs;
  }
}

Double_t CompSmearing(Double_t inputComponent){
  // Smears the components Ns and Nf according to fComponentRes
  //

  gRandom = new TRandom3();
  gRandom->SetSeed(1);

  if (fComponentRes == 0){
    return inputComponent;
  }else if (fComponentRes != 0 && inputComponent != 0){
    Double_t randomIs = gRandom->Gaus(0, inputComponent * fComponentRes / (sqrt(inputComponent)));
    return inputComponent + randomIs;
  }else{
    return inputComponent;
  }
}

void SwapVectorInt(vector <Int_t> *A, Int_t i, Int_t j){

  Int_t temp;

  temp = A->at(i);
  A->at(i) = A->at(j);
  A->at(j) = temp;
  
}

void SwapVectorDouble(vector <Double_t> *A, Int_t i, Int_t j){

  Double_t temp;

  temp = A->at(i);
  A->at(i) = A->at(j);
  A->at(j) = temp;
  
}

void SortCALIFA(Int_t left, Int_t right, vector <Int_t> *CALIFA_Id, vector <Int_t> *CALIFA_PointsInId, vector <Double_t> *CALIFA_X, vector <Double_t> *CALIFA_Y, vector <Double_t> *CALIFA_Z, vector <Double_t> *CALIFA_Init_X, vector <Double_t> *CALIFA_Init_Y, vector <Double_t> *CALIFA_Init_Z, vector <Double_t> *CALIFA_Theta, vector <Double_t> *CALIFA_Phi, vector <Double_t> *CALIFA_Time, vector <Double_t> *CALIFA_Ek, vector <Double_t> *CALIFA_Eloss, vector <Double_t> *CALIFA_Init_Eloss, vector <Double_t> *CALIFA_Init_Nf, vector <Double_t> *CALIFA_Init_Ns, vector <Double_t> *CALIFA_Nf, vector <Double_t> *CALIFA_Ns){
  
  Int_t i = left, j = right;

  Double_t pivot = CALIFA_Ek->at((left + right) / 2);
  
  /* partition */
  while (i <= j) {
    while (CALIFA_Ek->at(i) > pivot)
      i++;
    while (CALIFA_Ek->at(j) < pivot)
      j--;
    if (i <= j) {

      SwapVectorInt(CALIFA_Id,i,j);
      SwapVectorInt(CALIFA_PointsInId,i,j);
      SwapVectorDouble(CALIFA_X,i,j);
      SwapVectorDouble(CALIFA_Y,i,j);
      SwapVectorDouble(CALIFA_Z,i,j);
      SwapVectorDouble(CALIFA_Init_X,i,j);
      SwapVectorDouble(CALIFA_Init_Y,i,j);
      SwapVectorDouble(CALIFA_Init_Z,i,j);
      SwapVectorDouble(CALIFA_Theta,i,j);
      SwapVectorDouble(CALIFA_Phi,i,j);
      SwapVectorDouble(CALIFA_Time,i,j);
      SwapVectorDouble(CALIFA_Ek,i,j); 
      SwapVectorDouble(CALIFA_Eloss,i,j);
      SwapVectorDouble(CALIFA_Init_Eloss,i,j);
      SwapVectorDouble(CALIFA_Init_Nf,i,j);
      SwapVectorDouble(CALIFA_Init_Ns,i,j);
      SwapVectorDouble(CALIFA_Nf,i,j);
      SwapVectorDouble(CALIFA_Ns,i,j);


      i++;
      j--;
    }
  };
  /* recursion */
  if (left < j)
    SortCALIFA(left, j, CALIFA_Id, CALIFA_PointsInId, CALIFA_X, CALIFA_Y, CALIFA_Z, CALIFA_Init_X, CALIFA_Init_Y, CALIFA_Init_Z, CALIFA_Theta, CALIFA_Phi, CALIFA_Time, CALIFA_Ek, CALIFA_Eloss, CALIFA_Init_Eloss, CALIFA_Init_Nf, CALIFA_Init_Ns, CALIFA_Nf, CALIFA_Ns);
  if (i < right)
    SortCALIFA(i, right, CALIFA_Id, CALIFA_PointsInId, CALIFA_X, CALIFA_Y, CALIFA_Z, CALIFA_Init_X, CALIFA_Init_Y, CALIFA_Init_Z, CALIFA_Theta, CALIFA_Phi, CALIFA_Time, CALIFA_Ek, CALIFA_Eloss, CALIFA_Init_Eloss, CALIFA_Init_Nf, CALIFA_Init_Ns, CALIFA_Nf, CALIFA_Ns);
}


void SortCALIFA_OLD(vector <Int_t> *CALIFA_Id, vector <Int_t> *CALIFA_PointsInId, vector <Double_t> *CALIFA_X, vector <Double_t> *CALIFA_Y, vector <Double_t> *CALIFA_Z, vector <Double_t> *CALIFA_Init_X, vector <Double_t> *CALIFA_Init_Y, vector <Double_t> *CALIFA_Init_Z, vector <Double_t> *CALIFA_Theta, vector <Double_t> *CALIFA_Phi, vector <Double_t> *CALIFA_Time, vector <Double_t> *CALIFA_Ek, vector <Double_t> *CALIFA_Eloss, vector <Double_t> *CALIFA_Init_Eloss, vector <Double_t> *CALIFA_Init_Nf, vector <Double_t> *CALIFA_Init_Ns, vector <Double_t> *CALIFA_Nf, vector <Double_t> *CALIFA_Ns){

  // for(Int_t i=0 ; i<CALIFA_Id->size() ; i++){
  //   cout << CALIFA_Ek->at(i) << endl;
  // }
  
  for(Int_t i=0 ; i<CALIFA_Id->size()-1 ; i++){
    for(Int_t j=i+1 ; j<CALIFA_Id->size() ; j++){

      if(CALIFA_Ek->at(i)<CALIFA_Ek->at(j)){
	
	SwapVectorInt(CALIFA_Id,i,j);
	SwapVectorInt(CALIFA_PointsInId,i,j);
	SwapVectorDouble(CALIFA_X,i,j);
	SwapVectorDouble(CALIFA_Y,i,j);
	SwapVectorDouble(CALIFA_Z,i,j);
	SwapVectorDouble(CALIFA_Init_X,i,j);
	SwapVectorDouble(CALIFA_Init_Y,i,j);
	SwapVectorDouble(CALIFA_Init_Z,i,j);
	SwapVectorDouble(CALIFA_Theta,i,j);
	SwapVectorDouble(CALIFA_Phi,i,j);
	SwapVectorDouble(CALIFA_Time,i,j);
	SwapVectorDouble(CALIFA_Ek,i,j); 
	SwapVectorDouble(CALIFA_Eloss,i,j);
	SwapVectorDouble(CALIFA_Init_Eloss,i,j);
	SwapVectorDouble(CALIFA_Init_Nf,i,j);
	SwapVectorDouble(CALIFA_Init_Ns,i,j);
	SwapVectorDouble(CALIFA_Nf,i,j);
	SwapVectorDouble(CALIFA_Ns,i,j);

      }

    }
  }
  // cout << "Ordered " << endl;
  // for(Int_t i=0 ; i<CALIFA_Id->size() ; i++){
  //   cout << CALIFA_Ek->at(i) << endl;
  // }
}

void EraseVectorInt(vector <Int_t> *A, Int_t i){
  
  A->erase(A->begin()+i);
  
}

void EraseVectorDouble(vector <Double_t> *A, Int_t i){
  
  A->erase(A->begin()+i);
  
}

Double_t GetClosestApproachPoint_Bis(TVector3 Det_a0, TVector3 Det_a1, TVector3 Det_b0, TVector3 Det_b1, Double_t &Vertex_X, Double_t &Vertex_Y, Double_t &Vertex_Z){

  Double_t AX1 = (Det_a1.X()-Det_a0.X())/(Det_a1.Z()-Det_a0.Z());
  Double_t BX1 = Det_a1.X()-AX1*Det_a1.Z();
  Double_t AX2 = (Det_b1.X()-Det_b0.X())/(Det_b1.Z()-Det_b0.Z());
  Double_t BX2 = Det_b1.X()-AX2*Det_b1.Z();
  Double_t AY1 = (Det_a1.Y()-Det_a0.Y())/(Det_a1.Z()-Det_a0.Z());
  Double_t BY1 = Det_a1.Y()-AY1*Det_a1.Z();
  Double_t AY2 = (Det_b1.Y()-Det_b0.Y())/(Det_b1.Z()-Det_b0.Z());
  Double_t BY2 = Det_b1.Y()-AY2*Det_b1.Z();

  // double a1 = p[0];
  // double a2 = p[2];
  // double b1 = p[1];
  // double b2 = p[3];
  // double ap1 = pp[0];
  // double ap2 = pp[2];
  // double bp1 = pp[1];
  // double bp2 = pp[3];

  double a1 = BX1;
  double a2 = BY1;
  double b1 = AX1;
  double b2 = AY1;
  double ap1 = BX2;
  double ap2 = BY2;
  double bp1 = AX2;
  double bp2 = AY2;

  
  //cout << "TestFF :" << a1 << "," << a2 << "," << b1 << endl;
  //cout << "TestFFerr :" << ap1 << "," << ap2 << "," << bp1 << endl;

  double alpha, beta, A, B, C;
    
  alpha = (bp1*(a1-ap1)+bp2*(a2-ap2))/(bp1*bp1 + bp2*bp2 + 1);
  beta = (bp1*b1+bp2*b2+1)/(bp1*bp1 + bp2*bp2 + 1);
    
  A = beta*(bp1*bp1 + bp2*bp2 + 1) - (bp1*b1 + bp2*b2 + 1);
  B = (b1*b1 + b2*b2 + 1) - beta*(bp1*b1+bp2*b2+1);
  C = beta*(bp1*(ap1-a1) + bp2*(ap2-a2)) - (b1*(ap1-a1) + b2*(ap2-a2));
    

  //cout << "TestFF2 :" << A << "," << B << "," << C << endl;
  double sol1, solf1;
  double x,y,z,xp,yp,zp;
    
    
  sol1 = -(A*alpha + C)/(A*beta + B);
  solf1 = alpha + beta* sol1;

  //cout << "TestFF3 :" << sol1 << "," << solf1 <<  endl;
    
  x = a1 + b1*sol1;
  y = a2 + b2*sol1;
  z = sol1;
  xp = ap1 + bp1*solf1;
  yp = ap2 + bp2*solf1;
  zp = solf1;
    
  // xv = (x+xp)/2.;
  // yv = (y+yp)/2.;
  // zv = (z+zp)/2.;

  Vertex_X = (x+xp)/2.;
  Vertex_Y = (y+yp)/2.;
  Vertex_Z = (z+zp)/2.;

  return sqrt(pow((x-xp),2) + pow((y-yp),2) + pow((z-zp),2));
  
  /* cout << "Vertex 1st :" << x << "," << y << "," << z << endl;
     cout << "Vertex 2nd :" << xp << "," << yp << "," << zp << endl;
     cout << "Vertex middle :" << xv << "," << yv << "," << zv << endl;
    
     cout << "min dist " << sqrt(pow((x-xp),2) + pow((y-yp),2) + pow((z-zp),2)) << endl;*/
  
}

Double_t GetClosestApproachPoint(TVector3 a0, TVector3 a1, TVector3 b0, TVector3 b1, Double_t &Vertex_X, Double_t &Vertex_Y, Double_t &Vertex_Z){

   // Calculate denomitator
   TVector3 A = a1 - a0;
   TVector3 B = b1 - b0;
   Double_t magA = A.Mag();
   Double_t magB = B.Mag();

   Double_t scale_A = 1./magA;
   Double_t scale_B = 1./magB;
   TVector3 _A = scale_A*A;
   TVector3 _B = scale_B*B;

   TVector3 cross = _A.Cross(_B);
   Double_t denom = cross.Mag()*cross.Mag();


   // If lines are parallel (denom=0) test if lines overlap.
   // If they don't overlap then there is a closest point solution.
   // If they do overlap, there are infinite closest positions, but there is a closest distance
   if (denom == 0){
     cout << "Parallel trajectories" << endl;
     return 9999.;
   }

   // Lines criss-cross: Calculate the projected closest points
   TVector3 t = (b0 - a0);
   Double_t detA = t.X()*(_B.Y()*cross.Z()-_B.Z()*cross.Y()) - _B.X()*(t.Y()*cross.Z()-t.Z()*cross.Y()) + cross.X()*(t.Y()*_B.Z()-t.Z()*_B.Y()); //Determinant(t, _B, cross);
   Double_t detB = t.X()*(_A.Y()*cross.Z()-_A.Z()*cross.Y()) - _A.X()*(t.Y()*cross.Z()-t.Z()*cross.Y()) + cross.X()*(t.Y()*_A.Z()-t.Z()*_A.Y()); //Determinant(t, _A, cross);

   Double_t t0 = detA/denom;
   Double_t t1 = detB/denom;

   TVector3 pA = a0 + (_A * t0); // Projected closest point on segment A
   TVector3 pB = b0 + (_B * t1); // Projected closest point on segment B

   TVector3 ClosestPointsDist = pA-pB;
   Double_t result = ClosestPointsDist.Mag();

   Vertex_X = 0.5*(pA.X()+pB.X());
   Vertex_Y = 0.5*(pA.Y()+pB.Y());
   Vertex_Z = 0.5*(pA.Z()+pB.Z());
   
   return result;
 }

void run_proton_tracking(Char_t* InputFile, Bool_t MakeCalMap){

  Double_t FOOT_POS_SIG = 0.01500; //0.0150; //150um = 0.0150cm
  
  cout << "\n\n\t*************************************************" << endl;
  cout << "\t*                                               *" << endl;
  cout << "\t*       Proton tracking macros in R3BRoot       *" << endl;
  cout << "\t*                                               *" << endl;
  cout << "\t*************************************************" << endl;
  
  if(MakeCalMap){
    cout << ":::::::::::: Generating CALIFA Map ::::::::::::" << endl;
  }
  
  gRandom = new TRandom3();
  gRandom->SetSeed(1);

  //CALIFA ID_TO_POSITION Variables
  Double_t CALIFA_IDtoX[2432];
  Double_t CALIFA_IDtoY[2432];
  Double_t CALIFA_IDtoZ[2432];
  
  //Load CALIFA Map
  ifstream infile("CalFile/CALIFA_Mapping.dat");
  Int_t a;
  Double_t b, c, d;

  while (infile >> a >> b >> c >> d){
    //cout << a << "   " << b << "   " << c << "   " << d << endl; 
    CALIFA_IDtoX[a-1] = b;
    CALIFA_IDtoY[a-1] = c;
    CALIFA_IDtoZ[a-1] = d;
    //cout << a << "   " << CALIFA_IDtoX[a-1] << "   " << CALIFA_IDtoY[a-1] << "   " << CALIFA_IDtoZ[a-1] << endl; 
  }

  infile.close();

  //Cut for Punch-Through events (Ns VS Nf 2D-plot)
  TFile *cutfile_PT = new TFile("CalFile/cut_PT.root");

  TCutG *Cut_PT;
  Cut_PT = (TCutG*)cutfile_PT->Get("Cut_PT");

  cutfile_PT->Close();
  
  TRandom3 *Rand = new TRandom3(0);
  TRandom3 *RandX = new TRandom3(0);
  TRandom3 *RandY = new TRandom3(0);
  TRandom3 *RandZ = new TRandom3(0);
  TRandom3 *RandXX = new TRandom3(0);
  TRandom3 *RandYY = new TRandom3(0);
  TRandom3 *RandZZ = new TRandom3(0);
  
  TString inputName = (TString)"input/"+InputFile;
  TFile *f = new TFile(inputName);
  TTree *tree = (TTree*)f->Get("evt");
  if(tree==NULL){ cout << "\nERROR! No specified tree in the file " << InputFile  << endl; return; }
  TString outputName = (TString)"output/"+InputFile;
  TFile* output = new TFile(outputName, "RECREATE");
  TTree *OutTree = new TTree("OutTree","OutTree");
  
  //--------- GetData from sim tree --------
  TClonesArray * TraPointArray = new TClonesArray("R3BTraPoint");
  TClonesArray * CrystalPointArray = new TClonesArray("R3BCalifaPoint");
  TClonesArray * MCTrackArray   = new TClonesArray("R3BMCTrack");
  tree->SetBranchAddress("CrystalPoint",	&CrystalPointArray);
  tree->SetBranchAddress("TraPoint",	&TraPointArray);
  tree->SetBranchAddress("MCTrack",	&MCTrackArray);

  cout << "\n-- Input file: " << InputFile;
  cout << "\n-- Number of events in the file: " << tree->GetEntries() << endl;

  const Int_t Nentries=tree->GetEntries();

  R3BMCTrack *MCTrack;
  R3BTraPoint *TraPoint;
  R3BCalifaPoint *CrystalPoint;
  TVector3 FOOT_Pos;
  Double_t Vertex_X;
  Double_t Vertex_Y;
  Double_t Vertex_Z;
  vector <Int_t> *FOOT_Id = new vector <Int_t>;
  vector <Double_t> *FOOT_X = new vector <Double_t>;
  vector <Double_t> *FOOT_Y = new vector <Double_t>;
  vector <Double_t> *FOOT_Z = new vector <Double_t>;
  vector <Double_t> *FOOT_Theta = new vector <Double_t>;
  vector <Double_t> *FOOT_Phi = new vector <Double_t>;
  Int_t FOOT_Mult;
  Double_t ClosestDist;
  Bool_t FOOT2p = false;
  Bool_t FOOT1p = false;
  Double_t p1_Theta;
  Double_t p2_Theta;
  Double_t p1_Phi;
  Double_t p2_Phi;
  Double_t Angle_Opening;
  TVector3 Vertex;

  OutTree->Branch("FOOT_Id",&FOOT_Id);
  OutTree->Branch("FOOT_X",&FOOT_X);
  OutTree->Branch("FOOT_Y",&FOOT_Y);
  OutTree->Branch("FOOT_Z",&FOOT_Z);
  OutTree->Branch("FOOT_Theta",&FOOT_Theta);
  OutTree->Branch("FOOT_Phi",&FOOT_Phi);
  OutTree->Branch("FOOT_Mult",&FOOT_Mult);
  OutTree->Branch("ClosestDist",&ClosestDist);
  OutTree->Branch("Vertex_X",&Vertex_X);
  OutTree->Branch("Vertex_Y",&Vertex_Y);
  OutTree->Branch("Vertex_Z",&Vertex_Z);
  OutTree->Branch("p1_Theta",&p1_Theta);
  OutTree->Branch("p2_Theta",&p2_Theta);
  OutTree->Branch("p1_Phi",&p1_Phi);
  OutTree->Branch("p2_Phi",&p2_Phi);
  OutTree->Branch("Angle_Opening",&Angle_Opening);
  OutTree->Branch("FOOT2p",&FOOT2p);
  OutTree->Branch("FOOT1p",&FOOT1p);
  
  
  //CALIFA
  Int_t CALIFA_Mult;
  TVector3 CALIFA_Pos;
  vector <Int_t> *CALIFA_Id = new vector <Int_t>;
  vector <Int_t> *CALIFA_PointsInId = new vector <Int_t>;
  vector <Double_t> *CALIFA_X = new vector <Double_t>;
  vector <Double_t> *CALIFA_Y = new vector <Double_t>;
  vector <Double_t> *CALIFA_Z = new vector <Double_t>;
  vector <Double_t> *CALIFA_Init_X = new vector <Double_t>;
  vector <Double_t> *CALIFA_Init_Y = new vector <Double_t>;
  vector <Double_t> *CALIFA_Init_Z = new vector <Double_t>;
  vector <Double_t> *CALIFA_Theta = new vector <Double_t>;
  vector <Double_t> *CALIFA_Phi = new vector <Double_t>;
  vector <Double_t> *CALIFA_Time = new vector <Double_t>;
  vector <Double_t> *CALIFA_Ek = new vector <Double_t>; //After Punch-Through correction
  vector <Double_t> *CALIFA_Eloss = new vector <Double_t>; //Resolution && NONUniformity...
  vector <Double_t> *CALIFA_Init_Eloss = new vector <Double_t>; //Sim Output
  vector <Double_t> *CALIFA_Init_Nf = new vector <Double_t>; //fast light output [a.u]
  vector <Double_t> *CALIFA_Init_Ns = new vector <Double_t>; //slow light output [a.u]
  vector <Double_t> *CALIFA_Nf = new vector <Double_t>; //fast light output [a.u]
  vector <Double_t> *CALIFA_Ns = new vector <Double_t>; //slow light output [a.u]
  Double_t p1_Theta_CALIFA, p1_Phi_CALIFA;
  Double_t p2_Theta_CALIFA, p2_Phi_CALIFA;
  
  OutTree->Branch("CALIFA_Mult",&CALIFA_Mult);
  OutTree->Branch("CALIFA_Id",&CALIFA_Id);
  OutTree->Branch("CALIFA_PointsInId",&CALIFA_PointsInId);
  OutTree->Branch("CALIFA_X",&CALIFA_X);
  OutTree->Branch("CALIFA_Y",&CALIFA_Y);
  OutTree->Branch("CALIFA_Z",&CALIFA_Z);
  OutTree->Branch("CALIFA_Init_X",&CALIFA_Init_X);
  OutTree->Branch("CALIFA_Init_Y",&CALIFA_Init_Y);
  OutTree->Branch("CALIFA_Init_Z",&CALIFA_Init_Z);
  OutTree->Branch("CALIFA_Theta",&CALIFA_Theta);
  OutTree->Branch("CALIFA_Phi",&CALIFA_Phi);
  OutTree->Branch("CALIFA_Init_Eloss",&CALIFA_Init_Eloss);
  OutTree->Branch("CALIFA_Eloss",&CALIFA_Eloss);
  OutTree->Branch("CALIFA_Ek",&CALIFA_Ek);
  OutTree->Branch("CALIFA_Time",&CALIFA_Time);
  OutTree->Branch("CALIFA_Nf",&CALIFA_Nf);
  OutTree->Branch("CALIFA_Ns",&CALIFA_Ns);
  OutTree->Branch("CALIFA_Init_Nf",&CALIFA_Init_Nf);
  OutTree->Branch("CALIFA_Init_Ns",&CALIFA_Init_Ns);
  OutTree->Branch("p1_Theta_CALIFA",&p1_Theta_CALIFA);
  OutTree->Branch("p1_Phi_CALIFA",&p1_Phi_CALIFA);
  OutTree->Branch("p2_Theta_CALIFA",&p2_Theta_CALIFA);
  OutTree->Branch("p2_Phi_CALIFA",&p2_Phi_CALIFA);
  
  
  //MCTrack Info
  Double_t MCTrack_X, MCTrack_Y, MCTrack_Z;
  Int_t MCTrack_Mult;
  Double_t MCTrack1_Theta, MCTrack2_Theta;
  Double_t MCTrack1_Phi, MCTrack2_Phi;
  Double_t MCTrack1_E, MCTrack2_E;
  Double_t MCTrack1_Ek, MCTrack2_Ek;
  Double_t MCTrack1_Px, MCTrack2_Px;
  Double_t MCTrack1_Py, MCTrack2_Py;
  Double_t MCTrack1_Pz, MCTrack2_Pz;
  Double_t MCTrack_Angle_Opening;
  
  OutTree->Branch("MCTrack_X",&MCTrack_X);
  OutTree->Branch("MCTrack_Y",&MCTrack_Y);
  OutTree->Branch("MCTrack_Z",&MCTrack_Z);
  OutTree->Branch("MCTrack1_Theta",&MCTrack1_Theta);
  OutTree->Branch("MCTrack2_Theta",&MCTrack2_Theta);
  OutTree->Branch("MCTrack1_Phi",&MCTrack1_Phi);
  OutTree->Branch("MCTrack2_Phi",&MCTrack2_Phi);
  OutTree->Branch("MCTrack_Angle_Opening",&MCTrack_Angle_Opening);
  OutTree->Branch("MCTrack1_E",&MCTrack1_E);
  OutTree->Branch("MCTrack2_E",&MCTrack2_E);
  OutTree->Branch("MCTrack1_Ek",&MCTrack1_Ek);
  OutTree->Branch("MCTrack2_Ek",&MCTrack2_Ek);
  OutTree->Branch("MCTrack1_Px",&MCTrack1_Px);
  OutTree->Branch("MCTrack2_Px",&MCTrack2_Px);
  OutTree->Branch("MCTrack1_Py",&MCTrack1_Py);
  OutTree->Branch("MCTrack2_Py",&MCTrack2_Py);
  OutTree->Branch("MCTrack1_Pz",&MCTrack1_Pz);
  OutTree->Branch("MCTrack2_Pz",&MCTrack2_Pz);

  //Needed variables
  TVector3 p1, p2;
  
  TH1D *h1_CALIFA_X[2500];
  TH1D *h1_CALIFA_Y[2500];
  TH1D *h1_CALIFA_Z[2500];
  
  if(MakeCalMap){ //Making CALIFA Mapping

    //create histograms for each crystal
    for(Int_t i=0 ; i<2500 ; i++){
      char name[100];
      int f;
      f = sprintf(name,"h1_CALIFA_X_%i",i);
      h1_CALIFA_X[i] = new TH1D(name,name,1000,-60,60);
      f = sprintf(name,"h1_CALIFA_Y_%i",i);
      h1_CALIFA_Y[i] = new TH1D(name,name,1000,-60,60);
      f = sprintf(name,"h1_CALIFA_Z_%i",i);
      h1_CALIFA_Z[i] = new TH1D(name,name,1000,-80,80);
    }
  }

  
  cout << "\n-- Collecting track data.....\n";
  Int_t counter = 0;
  for(Int_t ev=0; ev<Nentries; ev++){
    if (ev%100==0) cout << "\r-- Working on entry " << ev << flush;

    TraPointArray->Clear();
    CrystalPointArray->Clear();
    MCTrackArray->Clear();

    FOOT_Id->clear();
    FOOT_X->clear();
    FOOT_Y->clear();
    FOOT_Z->clear();
    FOOT_Theta->clear();
    FOOT_Phi->clear();
    
    CALIFA_Id->clear();
    CALIFA_PointsInId->clear();
    CALIFA_X->clear();
    CALIFA_Y->clear();
    CALIFA_Z->clear();
    CALIFA_Init_X->clear();
    CALIFA_Init_Y->clear();
    CALIFA_Init_Z->clear();
    CALIFA_Theta->clear();
    CALIFA_Phi->clear();
    CALIFA_Init_Eloss->clear();
    CALIFA_Eloss->clear();
    CALIFA_Ek->clear();
    CALIFA_Time->clear();
    CALIFA_Nf->clear();
    CALIFA_Ns->clear();
    CALIFA_Init_Nf->clear();
    CALIFA_Init_Ns->clear();

    p1_Theta = -9999.;
    p2_Theta = -9999.;
    p1_Phi = -9999.;
    p2_Phi = -9999.;
    Angle_Opening = -9999.;

    p1_Theta_CALIFA = -9999.;
    p1_Phi_CALIFA = -9999.;
    p2_Theta_CALIFA = -9999.;
    p2_Phi_CALIFA = -9999.;
	
    MCTrack1_Theta = -9999.;
    MCTrack2_Theta = -9999.;
    MCTrack1_Phi = -9999.;
    MCTrack2_Phi = -9999.;
    MCTrack_Angle_Opening = -9999.;
    MCTrack1_E = -9999;
    MCTrack2_E = -9999;
    MCTrack1_Ek = -9999;
    MCTrack2_Ek = -9999;
    MCTrack1_Px = -9999;
    MCTrack2_Px = -9999;
    MCTrack1_Py = -9999;
    MCTrack2_Py = -9999;
    MCTrack1_Pz = -9999;
    MCTrack2_Pz = -9999;

    tree->GetEntry(ev);
      
    //MCTrack
    MCTrack_Mult = MCTrackArray->GetEntriesFast();
    for(Int_t i=0; i<MCTrackArray->GetEntriesFast(); i++){
      MCTrack = (R3BMCTrack*) MCTrackArray->At(i); 
      if(MCTrack->GetMotherId()==-1){
	MCTrack_X = MCTrack->GetStartX();
	MCTrack_Y = MCTrack->GetStartY();
	MCTrack_Z = MCTrack->GetStartZ();
	TVector3 Track_Mom;
	Track_Mom.SetXYZ(MCTrack->GetPx(),MCTrack->GetPy(),MCTrack->GetPz());
	if(i==0){
	  MCTrack1_E = MCTrack->GetEnergy();
	  MCTrack1_Px = MCTrack->GetPx();
	  MCTrack1_Py = MCTrack->GetPy();
	  MCTrack1_Pz = MCTrack->GetPz();
	  Double_t MCTrack_Gamma = MCTrack->GetEnergy()/MCTrack->GetMass();
	  MCTrack1_Ek = (MCTrack_Gamma-1.)*MCTrack->GetMass();
	  MCTrack1_Theta = Track_Mom.Theta();//*TMath::RadToDeg();
	  MCTrack1_Phi = Track_Mom.Phi();//*TMath::RadToDeg();
	}else{
	  MCTrack2_E = MCTrack->GetEnergy();
	  MCTrack2_Px = MCTrack->GetPx();
	  MCTrack2_Py = MCTrack->GetPy();
	  MCTrack2_Pz = MCTrack->GetPz();
	  Double_t MCTrack_Gamma = MCTrack->GetEnergy()/MCTrack->GetMass();
	  MCTrack2_Ek = (MCTrack_Gamma-1.)*MCTrack->GetMass();
	  MCTrack2_Theta = Track_Mom.Theta();//*TMath::RadToDeg();
	  MCTrack2_Phi = Track_Mom.Phi();//*TMath::RadToDeg();
	}
      }
      MCTrack_Angle_Opening = TMath::ACos(TMath::Sin(MCTrack1_Theta)*TMath::Sin(MCTrack2_Theta)*TMath::Cos(MCTrack2_Phi-MCTrack1_Phi)+TMath::Cos(MCTrack1_Theta)*TMath::Cos(MCTrack2_Theta));
    }

    //FOOT
    Bool_t IsDet1 = false;
    Bool_t IsDet2 = false;
    Bool_t IsDet3 = false;
    Bool_t IsDet4 = false;
    Bool_t IsDet5 = false;
    Bool_t IsDet6 = false;
    Bool_t IsDet7 = false;
    Bool_t IsDet8 = false;
    Bool_t IsDet9 = false;
    Bool_t IsDet10 = false;
  
    TVector3 Det1_Pos, Det2_Pos, Det3_Pos, Det4_Pos, Det5_Pos, Det6_Pos, Det7_Pos, Det8_Pos, Det9_Pos, Det10_Pos;
       
    FOOT_Mult = TraPointArray->GetEntriesFast();
    for(Int_t i=0; i<TraPointArray->GetEntriesFast(); i++){
      TraPoint = (R3BTraPoint*) TraPointArray->At(i);
      TraPoint->PositionOut(FOOT_Pos);
      //Retrieving positions
      Double_t Det_Orientation = 0.;
      if(TraPoint->GetDetCopyID()==1 || TraPoint->GetDetCopyID()==2 || TraPoint->GetDetCopyID()==3 || TraPoint->GetDetCopyID()==4 || TraPoint->GetDetCopyID()==6 || TraPoint->GetDetCopyID()==7){
	Det_Orientation = 75.;
      }else if( TraPoint->GetDetCopyID()==8 || TraPoint->GetDetCopyID()==9){
	Det_Orientation = 75.;
      }
      Double_t PosRes = RandX->Gaus(0,FOOT_POS_SIG);
      Double_t X = TMath::Sin(TMath::DegToRad()*(90.-Det_Orientation))*(PosRes+FOOT_Pos.X()/TMath::Sin(TMath::DegToRad()*(90.-Det_Orientation)));
      Double_t Y = FOOT_Pos.Y()+RandY->Gaus(0,FOOT_POS_SIG);
      Double_t Z = FOOT_Pos.Z() + FOOT_Pos.X()/TMath::Tan(TMath::DegToRad()*(90.-Det_Orientation)) - X/TMath::Tan(TMath::DegToRad()*(90.-Det_Orientation));
      // Double_t X = FOOT_Pos.X()+RandX->Gaus(0,FOOT_POS_SIG);
      // Double_t Z = FOOT_Pos.Z()+RandZ->Gaus(0,FOOT_POS_SIG);
      FOOT_Pos.SetXYZ(X,Y,Z);
      FOOT_X->push_back(FOOT_Pos.X());
      FOOT_Y->push_back(FOOT_Pos.Y());
      FOOT_Z->push_back(FOOT_Pos.Z());
      FOOT_Id->push_back(TraPoint->GetDetCopyID());
      if(TraPoint->GetDetCopyID()==1 && !IsDet1){
	IsDet1=true;
	Det1_Pos.SetXYZ(X,Y,Z);
      }
      if(TraPoint->GetDetCopyID()==2 && !IsDet2){
	IsDet2=true;
	Det2_Pos.SetXYZ(X,Y,Z);
      }
      if(TraPoint->GetDetCopyID()==3 && !IsDet3){
	IsDet3=true;
	Det3_Pos.SetXYZ(X,Y,Z);
      }
      if(TraPoint->GetDetCopyID()==4 && !IsDet4){
	IsDet4=true;
	Det4_Pos.SetXYZ(X,Y,Z);
      }
      if(TraPoint->GetDetCopyID()==5 && !IsDet5){
	IsDet5=true;
	Det5_Pos.SetXYZ(X,Y,Z);
      }
      if(TraPoint->GetDetCopyID()==6 && !IsDet6){
	IsDet6=true;
	Det6_Pos.SetXYZ(X,Y,Z);
      }
      if(TraPoint->GetDetCopyID()==7 && !IsDet7){
	IsDet7=true;
	Det7_Pos.SetXYZ(X,Y,Z);
      }
      if(TraPoint->GetDetCopyID()==8 && !IsDet8){
	IsDet8=true;
	Det8_Pos.SetXYZ(X,Y,Z);
      }
      if(TraPoint->GetDetCopyID()==9 && !IsDet9){
	IsDet9=true;
	Det9_Pos.SetXYZ(X,Y,Z);
      }
      if(TraPoint->GetDetCopyID()==10 && !IsDet10){
	IsDet10=true;
	Det10_Pos.SetXYZ(X,Y,Z);
      }
    }

    FOOT2p = false;
    FOOT1p = false;
    p1.SetXYZ(-9999.,-9999.,-9999.);
    p2.SetXYZ(-9999.,-9999.,-9999.);

    if(IsDet1 && IsDet2 && IsDet3 && IsDet4 && IsDet6 && IsDet7 && IsDet8 && IsDet9){
      FOOT2p = true;
    }

    if(FOOT2p){

      TVector3 Det_Pos1, Det_Pos2, Det_Pos3, Det_Pos4;

      //XYXY
      // Double_t A1 = (Det9_Pos.Y()-Det6_Pos.Y())/(Det9_Pos.Z()-Det6_Pos.Z());
      // Double_t B1 = Det9_Pos.Y() - A1*Det9_Pos.Z();
      // Double_t A2 = (Det8_Pos.Y()-Det7_Pos.Y())/(Det8_Pos.Z()-Det7_Pos.Z());
      // Double_t B2 = Det8_Pos.Y() - A2*Det8_Pos.Z();
      
      // Det_Pos1.SetXYZ(Det1_Pos.X(),A1*Det1_Pos.Z()+B1,Det1_Pos.Z());
      // Det_Pos2.SetXYZ(Det2_Pos.X(),A2*Det2_Pos.Z()+B2,Det2_Pos.Z());
      // Det_Pos3.SetXYZ(Det3_Pos.X(),A2*Det3_Pos.Z()+B2,Det3_Pos.Z());
      // Det_Pos4.SetXYZ(Det4_Pos.X(),A1*Det4_Pos.Z()+B1,Det4_Pos.Z());

      //YXYX
      // Double_t A1 = (Det9_Pos.X()-Det6_Pos.X())/(Det9_Pos.Z()-Det6_Pos.Z());
      // Double_t B1 = Det9_Pos.X() - A1*Det9_Pos.Z();
      // Double_t A2 = (Det8_Pos.X()-Det7_Pos.X())/(Det8_Pos.Z()-Det7_Pos.Z());
      // Double_t B2 = Det8_Pos.X() - A2*Det8_Pos.Z();
      
      // Det_Pos1.SetXYZ(A1*Det1_Pos.Z()+B1,Det1_Pos.Y(),Det1_Pos.Z());
      // Det_Pos2.SetXYZ(A2*Det2_Pos.Z()+B2,Det2_Pos.Y(),Det2_Pos.Z());
      // Det_Pos3.SetXYZ(A2*Det3_Pos.Z()+B2,Det3_Pos.Y(),Det3_Pos.Z());
      // Det_Pos4.SetXYZ(A1*Det4_Pos.Z()+B1,Det4_Pos.Y(),Det4_Pos.Z());

      //YXXY
      Double_t A1 = (Det4_Pos.X()-Det6_Pos.X())/(Det4_Pos.Z()-Det6_Pos.Z());
      Double_t B1 = Det4_Pos.X() - A1*Det4_Pos.Z();
      Double_t A2 = (Det3_Pos.X()-Det7_Pos.X())/(Det3_Pos.Z()-Det7_Pos.Z());
      Double_t B2 = Det3_Pos.X() - A2*Det3_Pos.Z();
      
      Det_Pos1.SetXYZ(A1*Det1_Pos.Z()+B1,Det1_Pos.Y(),Det1_Pos.Z());
      Det_Pos2.SetXYZ(A2*Det2_Pos.Z()+B2,Det2_Pos.Y(),Det2_Pos.Z());
      Det_Pos3.SetXYZ(A2*Det8_Pos.Z()+B2,Det8_Pos.Y(),Det8_Pos.Z());
      Det_Pos4.SetXYZ(A1*Det9_Pos.Z()+B1,Det9_Pos.Y(),Det9_Pos.Z());
      
      // Det_Pos1.SetXYZ(Det1_Pos.X(),Det6_Pos.Y(),0.5*(Det1_Pos.Z()+Det6_Pos.Z()));
      // Det_Pos2.SetXYZ(Det2_Pos.X(),Det7_Pos.Y(),0.5*(Det2_Pos.Z()+Det7_Pos.Z()));
      // Det_Pos3.SetXYZ(Det3_Pos.X(),Det8_Pos.Y(),0.5*(Det3_Pos.Z()+Det8_Pos.Z()));
      // Det_Pos4.SetXYZ(Det4_Pos.X(),Det9_Pos.Y(),0.5*(Det4_Pos.Z()+Det9_Pos.Z()));
      // cout << "event" << endl;
      //ClosestDist = GetClosestApproachPoint(Det_Pos1,Det_Pos4,Det_Pos2,Det_Pos3,Vertex_X,Vertex_Y,Vertex_Z);
      // cout << Vertex_Z << endl;
      ClosestDist = GetClosestApproachPoint_Bis(Det_Pos1,Det_Pos4,Det_Pos2,Det_Pos3,Vertex_X,Vertex_Y,Vertex_Z);
      // cout << Vertex_Z << endl;
      Vertex.SetXYZ(Vertex_X,Vertex_Y,Vertex_Z);

      //Trajectories of protons and angles
      // p1 = Det_Pos4 - Det_Pos1;
      // p2 = Det_Pos3 - Det_Pos2;
      p1 = Det_Pos1 - Vertex;
      p2 = Det_Pos2 - Vertex;
      
      
      p1_Theta = p1.Theta();//*TMath::RadToDeg();
      p2_Theta = p2.Theta();//*TMath::RadToDeg();
      p1_Phi = p1.Phi();//*TMath::RadToDeg();
      p2_Phi = p2.Phi();//*TMath::RadToDeg();

      Angle_Opening = TMath::ACos(TMath::Sin(p1_Theta)*TMath::Sin(p2_Theta)*TMath::Cos(p2_Phi-p1_Phi)+TMath::Cos(p1_Theta)*TMath::Cos(p2_Theta));
            
    }else if(IsDet1 && IsDet4 && IsDet6 && IsDet9){

      FOOT1p = true;
      
      TVector3 BeamPoint1,BeamPoint2;
      BeamPoint1.SetXYZ(0.,0.,0.);
      BeamPoint2.SetXYZ(0.,0.,1.);

      TVector3 Det_Pos1, Det_Pos4;
      
      Double_t A1 = (Det9_Pos.Y()-Det6_Pos.Y())/(Det9_Pos.Z()-Det6_Pos.Z());
      Double_t B1 = Det9_Pos.Y() - A1*Det9_Pos.Z();
      
      Det_Pos1.SetXYZ(Det1_Pos.X(),A1*Det1_Pos.Z()+B1,Det1_Pos.Z());
      Det_Pos4.SetXYZ(Det4_Pos.X(),A1*Det4_Pos.Z()+B1,Det4_Pos.Z());

      // Det_Pos1.SetXYZ(Det1_Pos.X(),Det6_Pos.Y(),0.5*(Det1_Pos.Z()+Det6_Pos.Z()));
      // Det_Pos4.SetXYZ(Det4_Pos.X(),Det9_Pos.Y(),0.5*(Det4_Pos.Z()+Det9_Pos.Z()));
      
      ClosestDist = GetClosestApproachPoint(Det_Pos1,Det_Pos4,BeamPoint1,BeamPoint2,Vertex_X,Vertex_Y,Vertex_Z);
      
      p1 = Det_Pos4 - Det_Pos1;

      p1_Theta = p1.Theta();//*TMath::RadToDeg();
      p1_Phi = p1.Phi();//*TMath::RadToDeg();

      Vertex.SetXYZ(Vertex_X,Vertex_Y,Vertex_Z);
      
    }else if(IsDet2 && IsDet3 && IsDet7 && IsDet8){

      FOOT1p = true;
      
      TVector3 Det_Pos2, Det_Pos3;

      Double_t A2 = (Det8_Pos.Y()-Det7_Pos.Y())/(Det8_Pos.Z()-Det7_Pos.Z());
      Double_t B2 = Det8_Pos.Y() - A2*Det8_Pos.Z();
      
      Det_Pos2.SetXYZ(Det2_Pos.X(),A2*Det2_Pos.Z()+B2,Det2_Pos.Z());
      Det_Pos3.SetXYZ(Det3_Pos.X(),A2*Det3_Pos.Z()+B2,Det3_Pos.Z());
      
      // Det_Pos2.SetXYZ(Det2_Pos.X(),Det7_Pos.Y(),0.5*(Det2_Pos.Z()+Det7_Pos.Z()));
      // Det_Pos3.SetXYZ(Det3_Pos.X(),Det8_Pos.Y(),0.5*(Det3_Pos.Z()+Det8_Pos.Z()));
      
      TVector3 BeamPoint1,BeamPoint2;
      BeamPoint1.SetXYZ(0.,0.,0.);
      BeamPoint2.SetXYZ(0.,0.,1.);
      
      ClosestDist = GetClosestApproachPoint(Det_Pos2,Det_Pos3,BeamPoint1,BeamPoint2,Vertex_X,Vertex_Y,Vertex_Z);
      
      p1 = Det_Pos3 - Det_Pos2;

      p1_Theta = p1.Theta();//*TMath::RadToDeg();
      p1_Phi = p1.Phi();//*TMath::RadToDeg();

      Vertex.SetXYZ(Vertex_X,Vertex_Y,Vertex_Z);
      
    }else{
      Vertex.SetXYZ(0.,0.,0.); //Set to target center
    }

    for(Int_t i=0; i<FOOT_X->size() ; i++){
      TVector3 FOOT_Hit;
      FOOT_Hit.SetXYZ(FOOT_X->at(i),FOOT_Y->at(i),FOOT_Z->at(i));
      TVector3 Temp_Traj = FOOT_Hit - Vertex;
      FOOT_Theta->push_back(Temp_Traj.Theta()*TMath::RadToDeg());
      FOOT_Phi->push_back(Temp_Traj.Phi()*TMath::RadToDeg());
    }

    
    //CALIFA
    Bool_t p1InCALIFA = false;
    Bool_t p2InCALIFA = false;
    CALIFA_Mult = CrystalPointArray->GetEntriesFast();
    for(Int_t i=0; i<CrystalPointArray->GetEntriesFast(); i++){
      CrystalPoint = (R3BCalifaPoint*) CrystalPointArray->At(i);

      Int_t HitId = CrystalPoint->GetCrystalId();
      Bool_t PrimaryHit = true;
      Int_t PrimaryHitID = -1;
      
      for(Int_t j=0 ; j<CALIFA_Id->size() ; j++){
	if(HitId==CALIFA_Id->at(j)){
	  PrimaryHit = false;
	  PrimaryHitID = j;
	}
      }

      if(!PrimaryHit){ //Adding Eloss and Ns/Nf for Points in same Crystal
	CALIFA_Init_Eloss->at(PrimaryHitID) += CrystalPoint->GetEloss();
	CALIFA_Eloss->at(PrimaryHitID) += NUSmearing(CrystalPoint->GetEloss());
	CALIFA_Init_Nf->at(PrimaryHitID) += CrystalPoint->GetNf();
	CALIFA_Init_Ns->at(PrimaryHitID) += CrystalPoint->GetNs();
	CALIFA_Nf->at(PrimaryHitID) += CrystalPoint->GetNf();
	CALIFA_Ns->at(PrimaryHitID) += CrystalPoint->GetNs();
	CALIFA_PointsInId->at(PrimaryHitID) += 1; 
	if(CALIFA_Time->at(PrimaryHitID)>CrystalPoint->GetTime()){
	  CALIFA_Time->at(PrimaryHitID) = CrystalPoint->GetTime();
	}
      }
      
      if(PrimaryHit){

	CrystalPoint->PositionIn(CALIFA_Pos);

	CALIFA_Init_X->push_back(CALIFA_Pos.X());
	CALIFA_Init_Y->push_back(CALIFA_Pos.Y());
	CALIFA_Init_Z->push_back(CALIFA_Pos.Z());
	CALIFA_X->push_back(CALIFA_IDtoX[HitId]);
	CALIFA_Y->push_back(CALIFA_IDtoY[HitId]);
	CALIFA_Z->push_back(CALIFA_IDtoZ[HitId]);
      
	CALIFA_Id->push_back(HitId);
	CALIFA_PointsInId->push_back(1);
	CALIFA_Init_Eloss->push_back(CrystalPoint->GetEloss());
	CALIFA_Eloss->push_back(NUSmearing(CrystalPoint->GetEloss()));
	CALIFA_Time->push_back(CrystalPoint->GetTime());
	CALIFA_Init_Nf->push_back(CrystalPoint->GetNf());
	CALIFA_Init_Ns->push_back(CrystalPoint->GetNs());
	CALIFA_Nf->push_back(CrystalPoint->GetNf());
	CALIFA_Ns->push_back(CrystalPoint->GetNs());
	
	//TVector3 Temp_Traj = CALIFA_Pos - Vertex;
	TVector3 CALIFA_PosFromMap;
	CALIFA_PosFromMap.SetXYZ(CALIFA_IDtoX[HitId],CALIFA_IDtoY[HitId],CALIFA_IDtoZ[HitId]);
	
	//TVector3 Temp_Traj = CALIFA_Pos - Vertex;
	TVector3 Temp_Traj = CALIFA_Pos;

	Double_t Theta = Temp_Traj.Theta()*TMath::RadToDeg();
	Double_t Phi = Temp_Traj.Phi()*TMath::RadToDeg();
      
	CALIFA_Theta->push_back(Theta);
	CALIFA_Phi->push_back(Phi);
      
	if(!p1InCALIFA && TMath::Abs(p1_Theta-Theta)<1 && TMath::Abs(p1_Phi-Phi)<2){
	  p1_Theta_CALIFA = Theta;
	  p1_Phi_CALIFA = Phi;
	  p1InCALIFA = true;
	}

	if(!p2InCALIFA && TMath::Abs(p2_Theta-Theta)<1 && TMath::Abs(p2_Phi-Phi)<2){
	  p2_Theta_CALIFA = Theta;
	  p2_Phi_CALIFA = Phi;
	  p2InCALIFA = true;
	}
      }
    }
    CALIFA_Mult = CALIFA_Id->size();
    
    //Take Threshold into account
    for(Int_t i=0 ; i<CALIFA_Id->size() ; i++){

      if(CALIFA_Eloss->at(i)<fThreshold){

	EraseVectorInt(CALIFA_Id,i);
	EraseVectorInt(CALIFA_PointsInId,i);
	EraseVectorDouble(CALIFA_X,i);
	EraseVectorDouble(CALIFA_Y,i);
	EraseVectorDouble(CALIFA_Z,i);
	EraseVectorDouble(CALIFA_Init_X,i);
	EraseVectorDouble(CALIFA_Init_Y,i);
	EraseVectorDouble(CALIFA_Init_Z,i);
	EraseVectorDouble(CALIFA_Theta,i);
	EraseVectorDouble(CALIFA_Phi,i);
	EraseVectorDouble(CALIFA_Time,i);
	EraseVectorDouble(CALIFA_Eloss,i);
	EraseVectorDouble(CALIFA_Init_Eloss,i);
	EraseVectorDouble(CALIFA_Init_Nf,i);
	EraseVectorDouble(CALIFA_Init_Ns,i);
	EraseVectorDouble(CALIFA_Nf,i);
	EraseVectorDouble(CALIFA_Ns,i);

	i--; //Move the Loop Back one step
	
      }
    }
    CALIFA_Mult = CALIFA_Id->size();
   
    //Reconstruction of Energy for Punch-Through Events & Adding Resolution
    for(Int_t i=0 ; i<CALIFA_Id->size() ; i++){

      if(fResolution > 0){
	Double_t temp = CALIFA_Eloss->at(i);
	CALIFA_Eloss->at(i) = ExpResSmearing(temp);
      }
      
      if(fComponentRes > 0){
	Double_t temp = CALIFA_Nf->at(i);
	CALIFA_Nf->at(i) = CompSmearing(temp);

	temp = CALIFA_Ns->at(i);
	CALIFA_Ns->at(i) = CompSmearing(temp);
      }
      
      CALIFA_Ek->push_back(CALIFA_Eloss->at(i));
      
      if(Cut_PT->IsInside(CALIFA_Nf->at(i),CALIFA_Ns->at(i))){
	CALIFA_Ek->at(i) = ElossToEnergy_PunchThrough(CALIFA_Eloss->at(i));
      }
    }

    if(CALIFA_Id->size()!=CALIFA_PointsInId->size() || CALIFA_Id->size()!=CALIFA_X->size() || CALIFA_Id->size()!=CALIFA_Y->size() || CALIFA_Id->size()!=CALIFA_Z->size() || CALIFA_Id->size()!=CALIFA_Init_X->size() || CALIFA_Id->size()!=CALIFA_Init_Y->size() || CALIFA_Id->size()!=CALIFA_Init_Z->size() || CALIFA_Id->size()!=CALIFA_Theta->size() || CALIFA_Id->size()!=CALIFA_Phi->size() || CALIFA_Id->size()!=CALIFA_Time->size() || CALIFA_Id->size()!=CALIFA_Ek->size() || CALIFA_Id->size()!=CALIFA_Eloss->size() || CALIFA_Id->size()!=CALIFA_Init_Eloss->size() || CALIFA_Id->size()!=CALIFA_Init_Nf->size() || CALIFA_Id->size()!=CALIFA_Init_Ns->size() || CALIFA_Id->size()!=CALIFA_Nf->size() || CALIFA_Id->size()!=CALIFA_Ns->size()){
      cout << "Problem" << endl;
    }
    
    //Sorting CALIFA Events 
    if(CALIFA_Id->size()>1){
      //By Decreasing CALIFA_Ek
      //SortCALIFA(0,CALIFA_Id->size()-1,CALIFA_Id, CALIFA_PointsInId, CALIFA_X, CALIFA_Y, CALIFA_Z, CALIFA_Init_X, CALIFA_Init_Y, CALIFA_Init_Z, CALIFA_Theta, CALIFA_Phi, CALIFA_Time, CALIFA_Ek, CALIFA_Eloss, CALIFA_Init_Eloss, CALIFA_Init_Nf, CALIFA_Init_Ns, CALIFA_Nf, CALIFA_Ns);

      //By Decreasing CALIFA_Eloss
      SortCALIFA(0,CALIFA_Id->size()-1,CALIFA_Id, CALIFA_PointsInId, CALIFA_X, CALIFA_Y, CALIFA_Z, CALIFA_Init_X, CALIFA_Init_Y, CALIFA_Init_Z, CALIFA_Theta, CALIFA_Phi, CALIFA_Time, CALIFA_Eloss, CALIFA_Ek, CALIFA_Init_Eloss, CALIFA_Init_Nf, CALIFA_Init_Ns, CALIFA_Nf, CALIFA_Ns);
    }
    
    if(MakeCalMap){ //Making CALIFA Mapping      
      for(Int_t i=0 ; i<CALIFA_Id->size() ; i++){
	h1_CALIFA_X[CALIFA_Id->at(i)-1]->Fill(CALIFA_Init_X->at(i));
	h1_CALIFA_Y[CALIFA_Id->at(i)-1]->Fill(CALIFA_Init_Y->at(i));
	h1_CALIFA_Z[CALIFA_Id->at(i)-1]->Fill(CALIFA_Init_Z->at(i));
      }
    }//End of CALIFA Mapping 
    
    OutTree->Fill();
  }//End of Event Loop
  
  cout << " " << endl;

  //Writing CALIFA Mapping in File
  if(MakeCalMap){

    ofstream myfile;
    myfile.open("CalFile/CALIFA_Mapping_new.dat");

    for(Int_t i=0 ; i<2432 ; i++){
      myfile << i+1 << "   " << h1_CALIFA_X[i]->GetMean() << "   " << h1_CALIFA_Y[i]->GetMean() << "   " << h1_CALIFA_Z[i]->GetMean() << "\n";
    }
    myfile.close();
    cout << "::::::::::: Wrote New CALIFA MAP in 'CalFile/CALIFA_Mapping_new.dat' :::::::::::" << endl;
  }

  OutTree->Write();
  output->Close();

  //Delete

  delete FOOT_Id;
  delete FOOT_X;
  delete FOOT_Y;
  delete FOOT_Z;
  delete FOOT_Theta;
  delete FOOT_Phi;


  delete CALIFA_Id;
  delete CALIFA_PointsInId;
  delete CALIFA_X;
  delete CALIFA_Y;
  delete CALIFA_Z;
  delete CALIFA_Init_X;
  delete CALIFA_Init_Y;
  delete CALIFA_Init_Z;
  delete CALIFA_Theta;
  delete CALIFA_Phi;
  delete CALIFA_Time;
  delete CALIFA_Eloss;
  delete CALIFA_Init_Eloss;
  delete CALIFA_Ek;
  delete CALIFA_Init_Nf;
  delete CALIFA_Init_Ns;
  delete CALIFA_Nf;
  delete CALIFA_Ns;
}

//Main function
int main(Int_t argc, Char_t* argv[]){
  Char_t *in_filename=0;
  Bool_t NeedHelp = kTRUE;
  Bool_t MakeCalMap = kFALSE;
  
  if (argc > 1){
    for (Int_t i = 0; i < argc; i++){
      if (strncmp(argv[i],"--file=",7) == 0){
	in_filename = argv[i]+7;
	NeedHelp = kFALSE;
      }
      if (strncmp(argv[i],"--CalMap",8) == 0){
	MakeCalMap = kTRUE;
	NeedHelp = kFALSE;
      }
    }
  }
  if (NeedHelp){
    cout << "\nUse the following flags: ";
    cout << "\n--file=/path/to/your/file/filename.root   : (mandatory) simulation input file";
    cout << "  " << endl;
    return 0;
  }
  
  run_proton_tracking(in_filename,MakeCalMap);
  return 0;
}
