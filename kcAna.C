//OAmppcに対する、oldHodのpmtの時刻の遅れ(あるいは早まり)をみる
//	手順:mppc(ch問わず)と各々のpmtの全てのヒットについて、直前のヒットとの時間差pmt-mppcをヒストグラムに詰める
//		mppcに対してpmtが遅れているのか早まっているのかは未知である
//		だから、mppcにヒットがあれば直前のpmtとpmt-mppcをフィルし、pmtにヒットがあれば直前のmppcとpmt-mppcをフィルする
//		mppcでのあるヒットとpmtでのあるヒットが同一のイベントである場合、そのようなものを集めれば、その時間差に対応するビンのエントリーは多くなる
//		mppcでのあるヒットとpmtでのあるヒットが無関係な場合、その時間差は一様にばらつく
//		従って、mppcとpmtとの相対的な時間差に対応するビンに、明確なピークが立つはず

//差分を取るときに最初の1イベント目を飛ばしていないが、たかだか1イベントゆえにその影響は無視できる

#include <math.h>

#define xLow 0
#define yLow -20     //みたいclockcount差の最小値
#define yUp 50       //みたいclockcount差の最大値
#define chNumPmt 12  // bits size for SIG of PMT
#define chNumMppc 64 // bits size for SIG of MPPC
#define clock 5      // [ns]

void kcAna(string intputFIle)
{
  //--------ヒストグラムの定義--------
  TH2D *tdcDiffH; // MPPCからのTDCの差分
  Int_t xBins = chNumPmt;
  Int_t yBins = (Int_t)(abs(yLow) + yUp);
  tdcDiffH = new TH2D("tdcDiffH", "TdcDiffFromMppc", xBins, (Double_t)xLow, (Double_t)chNumPmt, yBins, (Double_t)yLow, (Double_t)yUp);
  tdcDiffH->SetTitle("TdcDiffFromMppc;PMT ch;TdcDiffFromMppc [CLK];Entries/Bin");

  TH1D *chListMppcH = new TH1D("chListMppcH", "ChListMppc", chNumMppc, 0., chNumMppc); // MPPC, 横軸ch, 縦軸total hits
  TH1D *chListPmtH = new TH1D("chListPmtH", "chListPmtH", chNumPmt, 0., chNumPmt);     // PMT, 横軸ch, 縦軸total hits

  TH2D *hitMapMppcH = new TH2D("hitMapMppcH", "HitMap NewHod.", 70, -35.5, 34.5, 6, -12, 12); // NewHod, hitmap
  hitMapMppcH->SetStats(0);

  TH1D *tdcP3MppcH[chNumMppc];     // TDC, MPPC, P3基準
  TH1D *tdcMRSyncMppcH[chNumMppc]; // TDC, MPPC, MRSync基準
  for (Int_t i = 0; i < chNumMppc; i++)
  {
    tdcP3MppcH[i] = new TH1D(Form("MppcCh %d P3", i), Form("MppcCh %d P3", i), 250, 0., 3.); // ~ 3 [s]
    tdcP3MppcH[i]->SetTitle(Form("MppcCh %d P3;TdcFromP3 [s];Entries/Bin", i));
    tdcMRSyncMppcH[i] = new TH1D(Form("MppcCh %d MRSync", i), Form("MppcCh %d MRSync", i), 250, 0., 6.); // ~ 6 [us]
    tdcMRSyncMppcH[i]->SetTitle(Form("MppcCh %d MRSync;TdcFromMRSync [us];Entries/Bin", i));
  }

  TH1D *tdcP3PmtH[chNumPmt];
  TH1D *tdcMRSyncPmtH[chNumPmt];
  for (Int_t i = 0; i < chNumPmt; i++)
  {
    tdcP3PmtH[i] = new TH1D(Form("PmtCh %d P3", i), Form("PmtCh %d P3", i), 250, 0., 3.); // ~ 3 [s]
    tdcP3PmtH[i]->SetTitle(Form("PmtCh %d P3;TdcFromP3 [s];Entries/Bin", i));
    tdcMRSyncPmtH[i] = new TH1D(Form("PmtCh %d MRSync", i), Form("PmtCh %d MRSync", i), 250, 0., 6.); // ~ 6 [us]
    tdcMRSyncPmtH[i]->SetTitle(Form("PmtCh %d MRSync;TdcFromMRSync [us];Entries/Bin", i));
  }

  TString outputPdfName;
  outputPdfName.Form("kcAna.pdf");

  TCanvas *c1 = new TCanvas(outputPdfName.Data(), outputPdfName.Data());
  c1->Print(outputPdfName + "[", "pdf");

  TFile *inputFile = new TFile(intputFIle.c_str());
  TTree *tree;
  tree = (TTree *)inputFile->Get("tree");

  //--------mppc BRANCH--------
  std::vector<Long_t> *mppc; // 64 bits
  mppc = 0;
  TBranch *bmppc;
  bmppc = 0;
  tree->SetBranchAddress("mppc", &mppc, &bmppc);

  //--------pmt BRANCH--------
  std::vector<Int_t> *pmt; // 12 bits
  pmt = 0;
  TBranch *bpmt;
  bpmt = 0;
  tree->SetBranchAddress("pmt", &pmt, &bpmt);

  //--------tdcCnt BRANCH--------
  std::vector<Int_t> *tdcCnt; // tdcCnt: TDC(MPPC or PMT) from MRSync
  tdcCnt = 0;
  TBranch *btdcCnt;
  btdcCnt = 0;
  tree->SetBranchAddress("tdcCnt", &tdcCnt, &btdcCnt);

  //--------tdcMrSync BRANCH--------
  Int_t tdcMrSync; // tdcMrSync: TDC(MRSync) from P3
  tree->SetBranchAddress("tdcMrSync", &tdcMrSync);

  Int_t tdcLatestPmt[chNumPmt] = {}; // index i means i-th bit from LSB in 12 bits (i = 0,1,...)
  Int_t tdcLatestMppc = 0;           // doesn't care about ch

  Int_t nEntry = tree->GetEntries();
  for (Int_t iEntry = 0; iEntry < nEntry; iEntry++)
  {

    tree->GetEntry(iEntry);
    //spillCntにiEntry行目の値をいれる

    Long64_t tEntry = tree->LoadTree(iEntry);
    //iEntry行目のエントリを取得
    bpmt->GetEntry(tEntry);
    bmppc->GetEntry(tEntry);
    btdcCnt->GetEntry(tEntry);
    //mppcの値のすべてを取得
    for (Int_t iComponent = 0; iComponent < (Int_t)tdcCnt->size(); iComponent++)
    {
      for (Long_t iChMppc = 0; iChMppc < chNumMppc; iChMppc++)
      {
        if ((mppc->at(iComponent) >> iChMppc) & 1)
        {
          tdcP3MppcH[iChMppc]->Fill((tdcCnt->at(iComponent) + tdcMrSync) * clock * pow(10, -9)); // [s]
          tdcMRSyncMppcH[iChMppc]->Fill((tdcCnt->at(iComponent)) * pow(10, -3));                 // [us]
        }
      }
    }

    /*
    for (Int_t iComponent = 0; iComponent < (Int_t)tdcCnt->size(); iComponent++)
    {

      for (Int_t iChPmt = 0; iChPmt < chNumPmt; iChPmt++)
      {
        // --------hit on pmt--------
        if ((pmt->at(iComponent) >> iChPmt) & 1)
        {
          tdcDiffH->Fill(iChPmt, (tdcCnt->at(iComponent) + tdcMrSync) - tdcLatestMppc);
          chListPmtH->Fill(iChPmt);
          tdcP3PmtH[iChPmt]->Fill((tdcCnt->at(iComponent) + tdcMrSync) * clock * pow(10, -9)); // [s]
          tdcMRSyncPmtH[iChPmt]->Fill((tdcCnt->at(iComponent)) * clock * pow(10, -6));         // [us]
        }
      }

      // --------hit on mppc--------
      if (mppc->at(iComponent) != 0)
      {

        for (Int_t iChPmt = 0; iChPmt < chNumPmt; ++iChPmt)
        {
          tdcDiffH->Fill(iChPmt, tdcLatestPmt[iChPmt] - (tdcCnt->at(iComponent) + tdcMrSync));
        }

        for (Int_t iChMppc = 0; iChMppc < chNumMppc; iChMppc++)
        {
          if ((mppc->at(iComponent) >> iChMppc) & 1)
          {
            chListMppcH->Fill(iChMppc);
          }
        }

        //mppcの値のすべてを取得
        for (UInt_t j = 0; j < mppc->size(); ++j)
        {
          if ((mppc->at(j) >> 46) & 1)
          //n chぶんを右シフトして捨てる
          {
            tdcP3MppcH[46]->Fill((tdcCnt->at(j) + tdcMrSync) * clock * pow(10, -9)); // [s]
            tdcMRSyncMppcH[46]->Fill((tdcCnt->at(j)) * pow(10, -6));                 // [us]
          }
        }
      }

      // 最新値を更新
      for (Int_t iChPmt = 0; iChPmt < chNumPmt; iChPmt++)
      {
        // --------hit on pmt--------
        if ((pmt->at(iComponent) >> iChPmt) & 1)
        {
          tdcLatestPmt[iChPmt] = tdcCnt->at(iComponent) + tdcMrSync;
        }
      }

      // 最新値を更新

      // --------hit on mppc--------
      if (mppc->at(iComponent) != 0)
      {
        tdcLatestMppc = tdcCnt->at(iComponent) + tdcMrSync;
      }
    }
    */
  }
  inputFile->Close();

  //-------- hitmap --------
  hitMapMppcH->Draw("colz");
  for (int i = 0; i < 5; i++)
  {
    TLine *lh;
    if (i == 1 || i == 3)
    {
      lh = new TLine(-28.5, -8 + 4 * i, 27.5, -8 + 4 * i);
    }
    else
    {
      lh = new TLine(-32.5, -8 + 4 * i, 31.5, -8 + 4 * i);
    }
    lh->SetLineWidth(1);
    lh->Draw();
  }
  for (int i = 0; i < 57; i++)
  {
    TLine *lh;
    if (i % 7 == 0)
    {
      lh = new TLine(-28.5 + i, -8, -28.5 + i, 8);
    }
    else
    {
      lh = new TLine(-28.5 + i, -4, -28.5 + i, 4);
    }
    lh->SetLineWidth(1);
    lh->Draw();
  }
  TLine *lh1 = new TLine(-32.5, -8, -32.5, 8);
  lh1->SetLineWidth(1);
  lh1->Draw();
  TLine *lh2 = new TLine(31.5, -8, 31.5, 8);
  lh2->SetLineWidth(1);
  lh2->Draw();
  c1->Print(outputPdfName, "pdf");

  gStyle->SetOptStat("ne");
  //n:ヒストグラム名、e:entry数を表示

  //--------colock差--------
  tdcDiffH->Draw("colz");
  c1->Print(outputPdfName, "pdf");

  c1->SetLogy();
  chListMppcH->Draw();
  c1->Print(outputPdfName, "pdf");
  chListPmtH->Draw();
  c1->Print(outputPdfName, "pdf");

  //--------個別のmppc--------
  for (int i = 0; i < chNumMppc; ++i)
  {
    tdcP3MppcH[i]->Draw();
    gStyle->SetOptStat("nemruo");
    //n:ヒストグラム名、e:entry数、m:平均値、e:RMS、u:アンダーフロー数、o:オーバーフロー数を表示
    c1->Print(outputPdfName, "pdf");
    tdcMRSyncMppcH[i]->Draw();
    gStyle->SetOptStat("nemruo");
    //n:ヒストグラム名、e:entry数、m:平均値、e:RMS、u:アンダーフロー数、o:オーバーフロー数を表示
    c1->Print(outputPdfName, "pdf");
  }

  //--------個別のpmt--------
  for (int i = 0; i < chNumPmt; ++i)
  {
    tdcP3PmtH[i]->Draw();
    gStyle->SetOptStat("nemruo");
    //n:ヒストグラム名、e:entry数、m:平均値、e:RMS、u:アンダーフロー数、o:オーバーフロー数を表示
    c1->Print(outputPdfName, "pdf");
    tdcMRSyncPmtH[i]->Draw();
    gStyle->SetOptStat("nemruo");
    //n:ヒストグラム名、e:entry数、m:平均値、e:RMS、u:アンダーフロー数、o:オーバーフロー数を表示
    c1->Print(outputPdfName, "pdf");
  }

  /* Int_t projecCh = 0; */
  /* TH1D *projecChH = tdcDiffH->ProjectionY("projecCh histo", projecCh + 1, projecCh + 1); */
  /* projecChH->Draw(); */

  c1->Print(outputPdfName + "]", "pdf");
}
