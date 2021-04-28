//OAmppcに対する、oldHodのpmtの時刻の遅れ(あるいは早まり)をみる
//	手順:mppc(ch問わず)と各々のpmtの全てのヒットについて、直前のヒットとの時間差pmt-mppcをヒストグラムに詰める
//		mppcに対してpmtが遅れているのか早まっているのかは未知である
//		だから、mppcにヒットがあれば直前のpmtとpmt-mppcをフィルし、pmtにヒットがあれば直前のmppcとpmt-mppcをフィルする
//		mppcでのあるヒットとpmtでのあるヒットが同一のイベントである場合、そのようなものを集めれば、その時間差に対応するビンのエントリーは多くなる
//		mppcでのあるヒットとpmtでのあるヒットが無関係な場合、その時間差は一様にばらつく
//		従って、mppcとpmtとの相対的な時間差に対応するビンに、明確なピークが立つはず

//差分を取るときに最初の1イベント目を飛ばしていないが、たかだか1イベントゆえにその影響は無視できる

#define xLow 0
#define yLow -20	//みたいclockcount差の最小値
#define yUp 50		//みたいclockcount差の最大値
#define pmtChNum 12 //pmtのchの総数
#define mppcChNum 64

void Ana_kc(string fin){
  //--------ヒストグラムの定義--------
  TH2D *timeDiffMapH;
  Int_t xBins = pmtChNum;
  Int_t yBins = (Int_t)(abs(yLow) + yUp);
  timeDiffMapH = new TH2D("timeDiffMapH", "timeDiffMap histo", xBins, (Double_t)xLow, (Double_t)pmtChNum, yBins, (Double_t)yLow, (Double_t)yUp);
  timeDiffMapH->SetTitle("timediff of pmt to mppc;pmt ch;timediff from mppc [clockcount];entries/bin");

  TH1D *mppc_chHist = new TH1D("mppc_chHist","ch (mppc)",65,-0.5,64.5);
  TH1D *pmt_chHist = new TH1D("pmt_chHist","ch (pmt)",16,-0.5,15.5);
  
  TH1D *mppc_TDC_P3[mppcChNum];
  TH1D *mppc_TDC_MR[mppcChNum];
  for(int i=0; i<mppcChNum; ++i){
    mppc_TDC_P3[i] = new TH1D(Form("mppc_ch%d_P3",i),Form("ch%d (mppc)",i),10,0,10);
    mppc_TDC_P3[i] -> SetXTitle("TDC from P3");
    mppc_TDC_P3[i] -> SetYTitle("count/10bin");
    mppc_TDC_MR[i] = new TH1D(Form("mppc_ch%d_MrSync",i),Form("ch%d (mppc)",i),10,0,10);
    mppc_TDC_MR[i] -> SetXTitle("TDC from MrSync");
    mppc_TDC_MR[i] -> SetYTitle("count/10bin");
  }
  TH1D *pmt_TDC_P3[pmtChNum];
  TH1D *pmt_TDC_MR[pmtChNum];
  for(int i=0; i<pmtChNum; ++i){
    pmt_TDC_P3[i] = new TH1D(Form("pmt_ch%d_P3",i),Form("ch%d (pmt)",i),10,0,10);
    pmt_TDC_P3[i] -> SetXTitle("TDC from P3");
    pmt_TDC_P3[i] -> SetYTitle("count/10bin");
    pmt_TDC_MR[i] = new TH1D(Form("pmt_ch%d_MrSync",i),Form("ch%d (pmt)",i),10,0,10);
    pmt_TDC_MR[i] -> SetXTitle("TDC from MrSync");
    pmt_TDC_MR[i] -> SetYTitle("count/10bin");
  }
  TString Name;
  Name.Form("kc705.pdf");
  TCanvas *c1 = new TCanvas(Name.Data(), Name.Data());
  c1 -> Print(Name+"[","pdf");

  TFile *inputFile = new TFile(fin.c_str());
  TTree *tree;
  tree = (TTree *)inputFile->Get("tree");

  //--------pmtブランチ--------
  std::vector<Int_t> *pmt;
  pmt = 0;
  TBranch *bpmt;
  bpmt = 0;
  tree->SetBranchAddress("pmt", &pmt, &bpmt);

  //--------mppcブランチ--------
  std::vector<Long_t> *mppc;
  //hpc(32 bit)|lpc(32 bit)
  mppc = 0;
  TBranch *bmppc;
  bmppc = 0;
  tree->SetBranchAddress("mppc", &mppc, &bmppc);

  //--------tdcCntブランチ--------
  std::vector<Int_t> *tdcCnt;
  //tdcCnt: tdcMrSyncを基準としたtdcのclockCount
  tdcCnt = 0;
  TBranch *btdcCnt;
  btdcCnt = 0;
  tree->SetBranchAddress("tdcCnt", &tdcCnt, &btdcCnt);

  //--------tdcMrSyncブランチ--------
  Int_t tdcMrSync;
  //tdcMrSync: P3が基準
  tree->SetBranchAddress("tdcMrSync", &tdcMrSync);

  Int_t pmtTdcLatest[pmtChNum] = {};
  //pmtの各chでの直近のtdcCnt, インデックスがそのままch番号に対応
  //すべての要素を0で初期化
  Int_t mppcTdcLatest = 0;
  //mppc(ch問わず)での直近のtdcCnt

  Int_t nEntry = tree->GetEntries();
  for (Int_t iEntry = 0; iEntry < nEntry; ++iEntry){
    tree->GetEntry(iEntry);

    Long64_t tentry = tree->LoadTree(iEntry);

    bpmt->GetEntry(tentry);
    bmppc->GetEntry(tentry);
    btdcCnt->GetEntry(tentry);

    for (Int_t iComponent = 0; iComponent < (Int_t)tdcCnt->size(); ++iComponent){
      //--------pmtについて直近のmppcとのclockcount差をフィル--------
      for (Int_t iPmtCh = 0; iPmtCh < pmtChNum; ++iPmtCh){
	if ((pmt->at(iComponent) >> iPmtCh) & 1){
	  pmt_chHist -> Fill(iPmtCh);
	  pmt_TDC_P3[iPmtCh] -> Fill(tdcCnt->at(iComponent));
	  pmt_TDC_MR[iPmtCh] -> Fill(tdcCnt->at(iComponent) + tdcMrSync);
	  timeDiffMapH->Fill(iPmtCh, (tdcCnt->at(iComponent) + tdcMrSync) - mppcTdcLatest);
	}
      }

      for (Int_t iMppcCh = 0; iMppcCh < mppcChNum; ++iMppcCh){
	if ((mppc->at(iComponent) >> iMppcCh) & 1){
	  mppc_chHist -> Fill(iMppcCh);
	  mppc_TDC_P3[iMppcCh] -> Fill(tdcCnt->at(iComponent));
	  mppc_TDC_MR[iMppcCh] -> Fill(tdcCnt->at(iComponent) + tdcMrSync);
	}
      }
      
      //--------mppcについて直近のpmtとのclockcount差をフィル--------
      if (mppc->at(iComponent) != 0){
	for (Int_t iPmtCh = 0; iPmtCh < pmtChNum; ++iPmtCh){
	  timeDiffMapH->Fill(iPmtCh, pmtTdcLatest[iPmtCh] - (tdcCnt->at(iComponent) + tdcMrSync));
	}
      }

      //--------pmtについてpmtTdcLatestを更新--------
      for (Int_t iPmtCh = 0; iPmtCh < pmtChNum; ++iPmtCh){
	if ((pmt->at(iComponent) >> iPmtCh) & 1){
	  pmtTdcLatest[iPmtCh] = tdcCnt->at(iComponent) + tdcMrSync;
	}
      }

      //--------mppcについてpmtTdcLatestを更新--------
      if (mppc->at(iComponent) != 0){
	mppcTdcLatest = tdcCnt->at(iComponent) + tdcMrSync;
      }
    }
  }
  inputFile->Close();


  //gPad->SetLogz();
  c1 -> SetLogz();
  gStyle->SetOptStat("ne");
  //n:ヒストグラム名、e:entry数を表示

  timeDiffMapH->Draw("colz");
  c1 -> Print(Name,"pdf");

  c1 -> SetLogy();
  mppc_chHist -> Draw();
  c1 -> Print(Name,"pdf");
  for(int i=0; i<mppcChNum; ++i){
    mppc_TDC_P3[i] -> Draw();
    c1 -> Print(Name,"pdf");
    mppc_TDC_MR[i] -> Draw();
    c1 -> Print(Name,"pdf");
  }

  pmt_chHist -> Draw();
  c1 -> Print(Name,"pdf");
  for(int i=0; i<pmtChNum; ++i){
    pmt_TDC_P3[i] -> Draw();
    c1 -> Print(Name,"pdf");
    pmt_TDC_MR[i] -> Draw();
    c1 -> Print(Name,"pdf");
  }

  /* Int_t projecCh = 0; */
  /* TH1D *projecChH = timeDiffMapH->ProjectionY("projecCh histo", projecCh + 1, projecCh + 1); */
  /* projecChH->Draw(); */
  
  c1 -> Print(Name+"]","pdf");
}
