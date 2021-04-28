#ifndef tool_hxx
#define tool_hxx

#include "TAxis.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
TString getFileName(TString value){
  TObjArray *arrDir = value.Tokenize("/");
  int ndir = arrDir->GetEntries();
  
  TString buf = arrDir->At(ndir-1)->GetName();
  TObjArray *arrCom = buf.Tokenize(".");
  int ncomp = arrCom->GetEntries();
  TString out = "";
  for(int icomp=0;icomp<ncomp-1;icomp++){
    out += arrCom->At(icomp)->GetName();
    if(icomp!=ncomp-2) out += ".";
  }
  std::cout << out << std::endl;
  return out;
}

void setHist(TH1D *h, int color=4){
  h->SetLineColor(color);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetLabelSize(0.04);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);

  h->GetXaxis()->SetNdivisions(505);
}

void setHist(TH2D *h, int color=4){
  h->SetLineColor(color);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetLabelSize(0.04);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);
  ((TGaxis*)h->GetZaxis())->SetMaxDigits(3);
  
  h->GetXaxis()->SetNdivisions(505);
}

#endif
