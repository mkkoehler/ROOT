
#include <iostream>
#include <string>
#include <fstream>

#include "TROOT.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TArrow.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TVector3.h"
#include "TComplex.h"
#include "TProfile.h"
#include"TLatex.h"

#include"TTree.h"

void ReadTree(){

  TString strCent = "20_40"; // "00_20" , "20_40"
  TString strFile ="tree_VISHNU_" + strCent  + ".root" ;
  
  TFile fIn(strFile);
  fIn.ls();

  TTree *tr = (TTree*)fIn.Get("tree");
  //  tr->Print();


  double Tmin = 0.155;
  double Tmax = 0.159;
  
  TH2D *hBetaR = new TH2D("hBetaR","hBetaR",300,0.,30.,100,0.,1.);
  TH2D *hBetaT = new TH2D("hBetaT","hBetaT",100,0.,0.5,100,0.,1.);
  
  tr->Draw("vT:TMath::Sqrt(x*x+y*y)>>hBetaR",Form("T>=%lf && T<=%lf && vT > 0.01",Tmin,Tmax));
  tr->Draw("vT:T>>hBetaT","TMath::Sqrt(x*x+y*y)>0.01");

  double Mbeta  = hBetaR->GetMean(2);// 2-> y axis
  double MbetaE = hBetaR->GetMeanError(2);// 2-> y axis

  TProfile *pBetaRx = hBetaR->ProfileX("pBetaRx", 0, 100); 
  TProfile *pBetaTx = hBetaT->ProfileX("pBetaTx", 2, 100);

  hBetaT->SetTitle(";T/GeV ; #beta_{T}"); 
  pBetaRx->SetMarkerStyle(8);
  pBetaTx->SetMarkerStyle(8);

  hBetaR->GetXaxis()->SetRangeUser(0.,10.);
  hBetaR->SetTitle(";r (fm) ; #beta_{T}");
  hBetaR->GetXaxis()->SetNdivisions(5);
  hBetaR->GetYaxis()->SetNdivisions(5);
  hBetaR->GetXaxis()->SetTitleSize(0.06);
  hBetaR->GetYaxis()->SetTitleSize(0.06);
  hBetaR->GetXaxis()->SetTitleOffset(0.9);
  hBetaR->GetYaxis()->SetTitleOffset(0.9); 
  hBetaR->GetXaxis()->SetLabelSize(0.05);
  hBetaR->GetYaxis()->SetLabelSize(0.05);
  

  
  TF1 *fFitR = new TF1("fFitR","[0] + [1]*TMath::Power(x,[2])",0.,8.);
  fFitR->FixParameter(0,0.);
  //  fFitR->FixParameter(1,1.);
  pBetaRx->Fit(fFitR,"0");

  fFitR->SetLineStyle(2);
  fFitR->SetLineWidth(4);

  
  TCanvas *cBetaT = new TCanvas("cBetaT","cBetaT",800,600);
  cBetaT->SetLogz();
  cBetaT->SetGrid(1,1);
  hBetaT->DrawCopy("colz");
  pBetaTx->DrawCopy("same");
  
  
  TCanvas *cBetaR = new TCanvas("cBetaR","cBetaR",800,600);
  cBetaR->SetLogz();
  //  cBetaR->SetGrid(1,1);
  hBetaR->DrawCopy("col");
  //  pBetaRx->DrawCopy("same");
  //  fFitR->Draw("same");

  
  TLatex latex;
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top
  latex.DrawLatex(.4,.9,"VISHNU");
  latex.DrawLatex(.4,.8,Form("Centrality %s  %",strCent));
  latex.DrawLatex(.4,.7,Form("%1.3lf < T < %1.3lf GeV",Tmin,Tmax));
  latex.DrawLatex(.4,.6,Form("#LT#beta_{T}#GT = %1.2lf #pm %1.2lf",Mbeta,MbetaE));
  latex.DrawLatex(.4,.5,Form("n = %1.2lf #pm %1.2lf",fFitR->GetParameter(2),fFitR->GetParError(2)));

  
  
  TCanvas *cBetaRFancy = new TCanvas("cBetaRFancy","cBetaRFancy",800,600);
  cBetaRFancy->SetLogz();
  //  cBetaR->SetGrid(1,1);
  hBetaR->DrawCopy("col");
  //  pBetaRx->DrawCopy("same");
  //  fFitR->Draw("same");

  
  TLatex latex2;
  latex2.SetTextSize(0.05);
  latex2.SetTextAlign(13);  //align at top
  latex2.DrawLatex(.4,.9,"VISHNU");
  latex2.DrawLatex(.4,.8,"Centrality 0-20 %");
  latex2.DrawLatex(.4,.7,"T = T_{c}^{SHM}");

 


  if(!gROOT->IsBatch())return;
  
  // create GIF
  float TminG   = 0. ;
  float TmaxG   = 0. ;
  float Tmax0G  = 0.4;
  //  int framesG   = 100;
  int framesG   = 100;
  float stepG   = Tmax0G/framesG;
  TH2D *h2dG = new TH2D("h2dG","h2dG",240,-1.2,1.2,240,-1.2,1.2);
  TH2D *h2dGRT = new TH2D("h2dGRT","h2dGRT",300,-30.,30.,300,-30.,30.);

  
  TCanvas *cTestG = new TCanvas("cTestG","cTestG",800,800);
  cTestG->SetGrid(1,1);
  cTestG->SetLogz();
  gSystem->Unlink("evolBetaT.gif"); // delete old file
  
  TCanvas *cTestG2 = new TCanvas("cTestG2","cTestG2",800,800);
  cTestG2->SetGrid(1,1);
  cTestG2->SetLogz();
  gSystem->Unlink("evolRT.gif"); // delete old file
  
  for(int i=0 ; i<framesG ; i++){
    TmaxG = Tmax0G - i*stepG;
    TminG = TmaxG - stepG;
    
    cout << "frame -> " << i << endl;
    cout << "TminG -> " << TminG << " Tmax -> " << TmaxG << endl;
    
    tr->Draw("vx:vy>>h2dG",Form("T > 0.001 && T>%lf && T < %lf", TminG,TmaxG));
    h2dG->SetTitle(Form("T=%1.3lf GeV ;vx;vy",TminG));
    h2dG->DrawCopy("col");

    cTestG->Modified();
    cTestG->Update();
    //    cTestG->Print("evol.gif+");
    //    cTestG->Print("evol.gif+");

    
      tr->Draw("x:y>>h2dGRT",Form("T > 0.001 && T>%lf && T < %lf", TminG,TmaxG));
      cTestG2->Modified();
      cTestG2->Update();
      
      h2dGRT->SetTitle(Form("T=%1.3lf GeV ;x;y",TminG));
      h2dGRT->DrawCopy("col");
      cTestG2->Print("evolRT.gif+");
    
    
  }//for loop


  
}// end of macro


/*

******************************************************************************
*Tree    :tree      : tree                                                   *
*Entries : 16125680 : Total =       905282268 bytes  File  Size =   33840398 *
*        :          : Tree compression factor =  26.91                       *
******************************************************************************
*Br    0 :x         : x/D                                                    *
*Entries : 16125680 : Total  Size=  129323406 bytes  File Size  =    1160746 *
*Baskets :     3610 : Basket Size=    3199488 bytes  Compression= 111.35     *
*............................................................................*
*Br    1 :y         : y/D                                                    *
*Entries : 16125680 : Total  Size=  129323406 bytes  File Size  =    4464416 *
*Baskets :     3610 : Basket Size=    3199488 bytes  Compression=  28.95     *
*............................................................................*
*Br    2 :vx        : vx/D                                                   *
*Entries : 16125680 : Total  Size=  129327020 bytes  File Size  =    5653116 *
*Baskets :     3610 : Basket Size=    3200000 bytes  Compression=  22.86     *
*............................................................................*
*Br    3 :vy        : vy/D                                                   *
*Entries : 16125680 : Total  Size=  129327020 bytes  File Size  =    5665540 *
*Baskets :     3610 : Basket Size=    3200000 bytes  Compression=  22.81     *
*............................................................................*
*Br    4 :vT        : vT/D                                                   *
*Entries : 16125680 : Total  Size=  129327020 bytes  File Size  =    5708554 *
*Baskets :     3610 : Basket Size=    3200000 bytes  Compression=  22.64     *
*............................................................................*
*Br    5 :rho       : rho/D                                                  *
*Entries : 16125680 : Total  Size=  129330634 bytes  File Size  =    5703982 *
*Baskets :     3610 : Basket Size=    3200000 bytes  Compression=  22.66     *
*............................................................................*
*Br    6 :T         : T/D                                                    *
*Entries : 16125680 : Total  Size=  129323406 bytes  File Size  =    5270603 *
*Baskets :     3610 : Basket Size=    3199488 bytes  Compression=  24.52     *
*............................................................................*

 */
