#include <arpa/inet.h>
#include <errno.h>
#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <cstdint>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "tool.hxx"

using namespace std;

#define N_MPPC_CH 64
#define N_PMT_CH  12
            
FILE *fp;
FILE *flog;
TFile *file;
TTree *tree;
unsigned char buf[13];
//unsigned char* buf;

// Header: 10.5, BoardId: 0.5: Spill count: 2  
// Em count: 2, Reserved: 0.5, Footer: 11
// Reserved: 0.5 Signal ch: 8.25, Reserved: 0.25, Counter: 4
int DATA_UNIT = 13;

int num_events = 0;
int spillbr=0;
map<int,int> get_em;

int           boardId;
int           spillCnt;
int           tdcMrSync;
int           tdcCnt;
int           tdcOffset;
int           mrSync;
int           pmt;
long          mppc;
int           emCnt; 
int           numOfHitMppc;
int           numOfHitPmt;
vector<int>  *vTdcCnt;
vector<int>  *vPmt;
vector<long> *vMppc;
vector<int>  *vNumOfHitMppc;
vector<int>  *vNumOfHitPmt;

int usage(void){
  string message = "Usage: ./binary2root data_file";
  cerr << message << endl;
  return 0;
}

void clear(){
  vTdcCnt->clear();
  vPmt->clear();
  vMppc->clear();
  vNumOfHitMppc->clear();
  vNumOfHitPmt->clear();
}

bool correctByte(int nbyte){
  if(nbyte==0){
    tree->SetBranchAddress("spillCnt",&spillbr);
    TBranch *embranch = tree->Branch("emCnt",     &emCnt,     "emCnt/I"    );
    Long64_t loop = tree->GetEntries();
    for(Long64_t ev=0; ev<loop; ++ev){
      tree->GetEntry(ev);
      emCnt = get_em[spillbr];
      embranch -> Fill();
    }
    tree->Write();
    file->Close();
    if(feof(fp)){
      cout << "here" << endl;
      return false;
    }else if(ferror(fp)){
      printf("\ntotal event = %d\nfinish to make root file\n",num_events-1);
      printf("fread\n\n");
      //err(1, "fread");
    }else{
      printf("\ntotal event = %d\nfinish to make root file\n",num_events-1);
      printf("unknown error\n\n");
      //errx(1, "unknown error");
    }
  }else if(nbyte != DATA_UNIT){
    TBranch *embranch = tree->Branch("emCnt",     &emCnt,     "emCnt/I"    );
    tree->SetBranchAddress("spillCnt",&spillbr);
    Long64_t loop = tree->GetEntries();
    for(Long64_t ev=0;ev<loop;++ev){
      tree->GetEntry(ev);
      emCnt = get_em[spillbr];
      embranch -> Fill();
    }
    tree->Write();
    file->Close();
    printf("\ntotal event = %d\nfinish to make root file\n",num_events);
    printf("short read: try to read %d but returns %d bytes\n\n", DATA_UNIT, nbyte);
    //errx(1, "short read: try to read %d but returns %d bytes", DATA_UNIT, nbyte);
    return false;
  }
  return true;
}

bool isHeader(){
  unsigned int *pHead = (unsigned int *)(buf);
  //unsigned long *pH = (unsigned long*)&buf[0];
  //printf("%08x %08x\n",(int)*pHead[0],(int)*pHead[1]);
  //printf("     %16lx\n",(long)*pH);
  if((int)*pHead == 0x67452301){
  //unsigned int Head = ntohs(*pHead);
  //if( Head == 0x01234567){
    //cout << "header" << endl;
    return true;
  }
  return false;
}

// bool isHeader(){
//   unsigned int *pHead[2] = {(unsigned int *)&buf[0], (unsigned int *)&buf[4]};
//   unsigned long *pH = (unsigned long*)&buf[0];
//   //printf("%08x %08x\n",(int)*pHead[0],(int)*pHead[1]);
//   //printf("     %16lx\n",(long)*pH);
//   if( ((int)*pHead[0]==0x1200b0ab) && ((int)*pHead[1]==0x12705634) ){
//     return true;
//   } 
//   return false;
// }

bool isFooter(){
  //cout << "Footer" << endl;
  unsigned int *pFoot = (unsigned int *)(buf);
  if( *pFoot==0xaaaaaaaa){
    //cout << "footer" << endl;
    //printf("footer : 0x");
    //for(int i=0;i<11;i++) printf("%x",buf[i+2]);
    //puts("");
    unsigned short *pEmCnt = (unsigned short *)(buf+6);
    emCnt = ntohs(*pEmCnt);
    get_em[spillCnt] = emCnt;
    //cout << dec << emCnt << endl;
    return true;
  }
  return false;
}

// bool isFooter(){
//   unsigned int   *pFoot[2]  = {(unsigned int   *)&buf[9], (unsigned int *)&buf[5]};
//   if( ((int)*pFoot[0]==0xAAAAAAAA) && ((int)*pFoot[1]==0xAAAAAAAA) ){
//     //printf("footer : 0x");                                                                                                                 
//     //for(int i=0;i<11;i++) printf("%x",buf[i+2]);                                                                                           
//     //puts("");                                                                                                                              
//     unsigned short *pEmCnt = (unsigned short *)&buf[0];
//     emCnt = ntohs(*pEmCnt);
//     return true;
//   }
//   return false;
// }

bool getMppcNHitCh(long val, int &num){
  num = 0;
  bool ok = true;
  for(int ich=0;ich<N_MPPC_CH;ich++){
    bool usedCh = (ich==32) | (ich==33) | (ich==42) | (ich==43) | (ich>=52);   // FIXME
    bool hit = (val >> ich) & 0x1;
    if(hit) num++;
    if(hit & (not usedCh) & ok) ok = false;
  } 
  return ok;
}

bool getPmtNHitCh(int val, int &num){
  num = 0;
  bool ok = true;
  for(int ich=0;ich<N_PMT_CH;ich++){
    bool hit = (val >> ich) & 0x1;
    if(hit) num++;
    if(hit & (ich==(N_PMT_CH-1))) ok = false;
  } 
  return ok;
}

int decode(){
  int nread = 0;
  while(1){
    int nbyte = fread(buf, 1, DATA_UNIT, fp);
    nread++;
    //for(int j=0; j<13; j++){
    //  cout << hex << +(unsigned char)buf[j] << " ";
    //}
    //cout << endl;
    if(not correctByte(nbyte)) {printf("incorrect header\n"); return -1;}
    //cout << nread << " " << nbyte << endl;
    if(not isHeader()) continue;
    
    // Read header
    unsigned char   *pBoardId   = (unsigned char  *)&buf[6];
    unsigned short  *pSpillCnt  = (unsigned short *)&buf[4];
    boardId = (int)*pBoardId & 0xF;
    spillCnt = ntohs(*pSpillCnt);

    tdcMrSync = 0;
    tdcOffset = 0;
    int preTdcCnt = -1;

    while(1){
      nbyte = fread(buf, 1, DATA_UNIT, fp);
      nread++;
      //cout << dec << nread << " ";
      // for(int j=0; j<13; j++){
      // 	cout << setw(2) <<  hex << +(unsigned char)buf[j] << " ";
      // }
      // cout << endl;
      //if (nread > 30000 && nread % 100 == 0) {
      //while (std::getchar() != '\n');
      //}
      //cout << nread << " " << nbyte << endl;
      // if(nread==384866 || nread==384865 || nread==384867){
      // 	for(int j=0; j<13; j++){
      // 	  cout << hex << +(unsigned char)buf[j] << " ";
      // 	}
      // 	cout << endl;

      if(not correctByte(nbyte)) {printf("incorrect event data or fotter\n"); return -1;}
      if(isFooter()) break;
      
      // Data format: {SIGNAL[63:0], OLDH[11:0], MR_SYNC, COUNTER[26:0]}
      unsigned int   *pMppc[2] = {(unsigned int   *)&buf[0], (unsigned int *)&buf[4]};
      unsigned short *pPmt     =  (unsigned short *)&buf[8];
      unsigned int   *pTdcCnt  =  (unsigned int   *)&buf[9];
      mppc   = ntohl(*pMppc[0]) * pow(2,4*8) + ntohl(*pMppc[1]);
      pmt    = ntohs(*pPmt) >> 4;
      mrSync = (ntohs(*pPmt) >> 3) & 0x1;
      tdcCnt = ntohl(*pTdcCnt) & 0x07FFFFFF;

      // Data checker
      int nMppcHit(0), nPmtHit(0);
      if( (not getMppcNHitCh(mppc,nMppcHit)) || (not getPmtNHitCh(pmt,nPmtHit)) || (tdcCnt>2048 && abs(tdcCnt+tdcOffset-preTdcCnt)>pow(2,16)) ){
        fprintf(flog," read : %8d  ",nread);
        for(int i=0;i<13;i++) fprintf(flog,"%02x",buf[i]);
        fprintf(flog,"\n");

        unsigned char *tmpBuf[2] = {(unsigned char *)&buf[11], (unsigned char *)&buf[10]};
        int nShift = 3;
        if(__builtin_popcount((int)*tmpBuf[0])>=3 )      nShift = 1;
        else if(__builtin_popcount((int)*tmpBuf[1])>=3 ) nShift = 2;
        else                                             nShift = 3;
        //fread(buf, 1, DATA_UNIT-nShift, fp);
        continue;
      }else if( (preTdcCnt<(pow(2,27)-2048) || tdcCnt>pow(2,12)) && (tdcCnt+tdcOffset)<preTdcCnt){
        fprintf(flog," read : %8d  ",nread);
        for(int i=0;i<13;i++) fprintf(flog,"%02x",buf[i]);
        fprintf(flog," %d - %d = %d\n",tdcCnt+tdcOffset,preTdcCnt,tdcCnt+tdcOffset-preTdcCnt);
        continue;
      }

      if((tdcCnt+tdcOffset)<preTdcCnt && preTdcCnt>(pow(2,27)-2048) && tdcCnt<pow(2,12)) tdcOffset+=pow(2,27);
      tdcCnt += tdcOffset;

      if(mrSync==1){
        if(vPmt->size()!=0) tree->Fill();
        clear();
        tdcMrSync = tdcCnt;
      }

      if(pmt!=0 || mppc!=0){
        vTdcCnt->push_back( tdcCnt - tdcMrSync );
        vPmt->push_back(pmt);
        vMppc->push_back(mppc);
        vNumOfHitPmt->push_back( nPmtHit );
        vNumOfHitMppc->push_back( nMppcHit );
      }
      //printf("%x, %x, %x\n",mrSync, coin, mppcCh);
      preTdcCnt = tdcCnt;
    }
    tree->Fill();
    num_events++;
    cout << "events: " << num_events << endl;
    clear();
  }
  tree->Fill();
  TBranch *embranch = tree->Branch("emCnt",     &emCnt,     "emCnt/I"    );
  tree->SetBranchAddress("spillCnt",&spillbr);
  Long64_t loop = tree->GetEntries();
  for(Long64_t ev=0;ev<loop;++ev){
    tree->GetEntry(ev);
    emCnt = get_em[spillbr];
    embranch -> Fill();
  }
  tree->Write();
  file->Close();
  return 0;
}

int main(int argc, char **argv){
  if(argc!=2){
    usage();
    exit(1);
  }

  TString input  = argv[1];
  TString output =  getFileName(input) + ".root";
  TString logfile = getFileName(input) + ".log";
  printf("Convert %s into %s (log file: %s)\n",input.Data(), output.Data(),logfile.Data());
  
  // Set read data size
  //buf = (unsigned char*)malloc(DATA_UNIT*sizeof(unsigned char));
  //if(buf==NULL){
  //  fprintf(stderr,"failed to malloc my_buf\n");
  //  exit(1);
  //}
  file = new TFile(output.Data(),"recreate");
  tree = new TTree("tree","kc705");

  vTdcCnt       = 0;
  vPmt          = 0;
  vMppc         = 0;
  vNumOfHitMppc = 0;
  vNumOfHitPmt  = 0;

  tree->Branch("boardId",   &boardId,   "boardId/I");
  tree->Branch("spillCnt",  &spillCnt,  "spillCnt/I" );
  tree->Branch("tdcMrSync", &tdcMrSync, "tdcMrSync/I");
  tree->Branch("tdcCnt",       &vTdcCnt       );
  tree->Branch("pmt",          &vPmt          );
  tree->Branch("mppc",         &vMppc         );
  tree->Branch("numOfHitMppc", &vNumOfHitMppc );
  tree->Branch("numOfHitPmt",  &vNumOfHitPmt  );
  
  vTdcCnt       = new vector<int>;
  vPmt          = new vector<int>;
  vMppc         = new vector<long>;
  vNumOfHitMppc = new vector<int>;
  vNumOfHitPmt  = new vector<int>;
  
  fp = fopen(input.Data(),"r");
  if(fp==NULL){
    err(1,"fopen (data file)");
  }

  flog = fopen(logfile.Data(),"w");
  if(flog==NULL){
    err(1,"fopen (log file)");
  }

  decode();
  printf("\ntotal event = %d\nfinish to make root file\n",num_events);
  return 0;
}
