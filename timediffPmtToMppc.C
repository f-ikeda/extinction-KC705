//mppcに対する、oldHodのpmtの時刻の遅れ(あるいは早まり)をみる
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

void timediffPmtToMppc()
{
	//--------canvasの定義--------
	TCanvas *Canvas = new TCanvas("c1", "canvas", 1);
	//("name","title",form=1デフォルトのサイズ)
	//あるいは("name","title",横,縦)
	Canvas->cd();

	//--------ヒストグラムの定義--------
	TH2F *timeDiffMapH;
	//x軸: pmtのch番号
	//y軸: mppcに対するpmtの、clockcount差
	Int_t xBins = pmtChNum;
	//x軸のbinの数
	Int_t yBins = (Int_t)(abs(yLow) + yUp);
	//y軸のbinの数
	//1 bin=1 clockcount
	timeDiffMapH = new TH2F("timeDiffMapH", "timeDiffMap histo", xBins, (Double_t)xLow, (Double_t)pmtChNum, yBins, (Double_t)yLow, (Double_t)yUp);
	timeDiffMapH->SetTitle("timediff of pmt to mppc;pmt ch;timediff from mppc [clockcount];entries/bin");

	//--------ファイルを開く--------
	TFile *inputFile;
	inputFile = TFile::Open("beamOn-10shoot-20201217.root");
	TTree *tree;
	tree = (TTree *)inputFile->Get("tree");

	//--------pmtブランチ--------
	std::vector<Int_t> *pmt;
	//pmt(12 bit)
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
	for (Int_t iEntry = 0; iEntry < nEntry; ++iEntry)
	{
		tree->GetEntry(iEntry);

		Long64_t tentry = tree->LoadTree(iEntry);

		bpmt->GetEntry(tentry);
		bmppc->GetEntry(tentry);
		btdcCnt->GetEntry(tentry);

		for (Int_t iComponent = 0; iComponent < (Int_t)tdcCnt->size(); ++iComponent)
		{

			//--------pmtについて直近のmppcとのclockcount差をフィル--------
			for (Int_t iPmtCh = 0; iPmtCh < pmtChNum; ++iPmtCh)
			//pmtの各chについて
			{
				if ((pmt->at(iComponent) >> iPmtCh) & 1)
				//iPmtChにヒットがあれば
				{
					timeDiffMapH->Fill(iPmtCh, (tdcCnt->at(iComponent) + tdcMrSync) - mppcTdcLatest);
					//iPmtChに、直近のmppcのヒットとのclockcount差をフィル
				}
			}

			//--------mppcについて直近のpmtとのclockcount差をフィル--------
			if (mppc->at(iComponent) != 0)
			//mppcのいずれかのchにヒットがあれば
			{
				for (Int_t iPmtCh = 0; iPmtCh < pmtChNum; ++iPmtCh)
				//pmtの各chについて
				{
					timeDiffMapH->Fill(iPmtCh, pmtTdcLatest[iPmtCh] - (tdcCnt->at(iComponent) + tdcMrSync));
					//iPmtChの直近のヒットとのclockcount差をフィル
				}
			}

			//--------pmtについてpmtTdcLatestを更新--------
			for (Int_t iPmtCh = 0; iPmtCh < pmtChNum; ++iPmtCh)
			//pmtの各chについて
			{
				if ((pmt->at(iComponent) >> iPmtCh) & 1)
				//iPmtChにヒットがあれば
				{
					pmtTdcLatest[iPmtCh] = tdcCnt->at(iComponent) + tdcMrSync;
				}
			}

			//--------mppcについてpmtTdcLatestを更新--------
			if (mppc->at(iComponent) != 0)
			//mppcのいずれかのchにヒットがあれば
			{
				mppcTdcLatest = tdcCnt->at(iComponent) + tdcMrSync;
			}
		}
	}

	gPad->SetLogz();
	//z軸をログスケールに

	gStyle->SetOptStat("ne");
	//n:ヒストグラム名、e:entry数を表示

	timeDiffMapH->Draw("colz");
	//z軸を色で表示
	//timeDiffMapH->Draw("lego2z");
	//z軸を立体的に色で表示

	/*何かがおかしい、動かない
	//--------X軸のあるchでスライスしてY軸に射影する--------
	Int_t projecCh = 0;
	//スライスしてみたいch
	TH1D *projecChH = timeDiffMapH->ProjectionY("projecCh histo", projecCh + 1, projecCh + 1);
	//ビン番号は1始まりであることに注意
	projecChH->Draw();
	*/

	inputFile->Close();
}
