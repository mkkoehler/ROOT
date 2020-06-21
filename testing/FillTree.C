
#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TArrow.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TVector3.h"
#include "TComplex.h"


#include"TTree.h"

void FillTree(){

  string cent = "20_40";// "00_20" , "20_40"
  string fn = "./VISHNUout/LHC_CENT";
  fn += cent ;
  fn += "_avg_Vx_Vy_T.txt";

  //  cout << "fn -> " << fn << endl;
  
  ifstream f(fn.c_str());


  TFile *fOut = new TFile(Form("tree_VISHNU_%s.root",cent.c_str()),"RECREATE");

  
  // read header:
  // number of time steps, number of cells in x and y
  // starting time, the maximum time,
  // the minimum and the maximum x and y values.

  const Double_t sqrts = 2760; // GeV
  const Double_t mp = 0.938;   // GeV
  // const Double_t Ymax = log(sqrts/mp);
  const Double_t Ymax = 8;
  Double_t vT, rho;
  Double_t vx, vy, T;
  Double_t x,y;
  Double_t tstart, tmax, xmin, xmax, ymin, ymax, dx2, dy2, dt2;
  Int_t nt, nx, ny;
  string dummy2;
  char line[255];

  f.getline(line,255);
  f >> dummy2 >> tstart >> xmin >> ymin ;
  f.getline(line,255);
  f.getline(line,255);
  f >> dummy2 >> tmax >> xmax >> ymax;
  f.getline(line,255);
  f.getline(line,255);
  f >> dummy2 >> dt2 >> dx2 >> dy2;
  f.getline(line,255);
  f.getline(line,255);
  f >> dummy2 >> nt >> nx >> ny;

  
  TTree *tree = new TTree("tree","tree");
  tree->Branch("x",&x);
  tree->Branch("y",&y);
  tree->Branch("vx",&vx);
  tree->Branch("vy",&vy);
  tree->Branch("vT",&vT);
  tree->Branch("rho",&rho);
  tree->Branch("T",&T);
  
  Double_t dx = (xmax-xmin)/(nx-1);
  Double_t dy = (ymax-ymin)/(ny-1);
  Double_t dt = (tmax-tstart)/(nt-1);

  // loop over time steps
  for (Int_t it=0; it <nt; it++) {

    // proper time tau
    Double_t t = tstart + it * dt;
    Double_t tfinal = tstart + nt * dt;

    // loop of x any y coordinates
    for (Int_t ix = 0; ix < nx; ix++) {

      x = xmin + ix * dx;

      for (Int_t iy = 0; iy < ny; iy++) {

        y = ymin + iy * dy;

        //Read in velocity in x and y, temperature T
        f >> vx >> vy >> T;
	
	// loop over space-time rapidity: dY
        const Int_t nYbins = 20;
        const Double_t dY = Ymax / nYbins;

        // only integrate for Y>0 (add factor 2 further below)
        for (Int_t iY = 0; iY < nYbins; iY++) {

          Double_t Y = iY * dY + dY/2.;

          // determine four-velocity vector of the fluid cell
          vT = sqrt(vx*vx + vy*vy);

          // transverse rapidity
          rho = 0.5 * log((1.+vT)/(1.-vT));

	  tree->Fill();
	  
	} //End loop Y
      } //End loop y
    } //End loop x
  } //End loop t


  tree->Write();
  fOut->Close();
  
  
}
