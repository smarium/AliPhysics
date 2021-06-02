/*An analysis task that plots the TPC dE/dx spectra Protons, Kaons, Pions and Antiprotons*/
/*N. Jacazio, S. Marium and A. Kalweit*/
/*01.06.2021*/

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TH1D.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TLatex.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliVVertex.h"
#include "AliExternalTrackParam.h"
#include "AliAODHeader.h"
#include "AliAODMCParticle.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisTaskTpcSpectra.h"
#include "AliMultSelection.h"

class AliAnalysisTaskTpcSpectra;

using namespace std;

ClassImp(AliAnalysisTaskTpcSpectra)

  //____________________________________________________________________________________//

  AliAnalysisTaskTpcSpectra::AliAnalysisTaskTpcSpectra() : AliAnalysisTaskTpcSpectra("stdtask")
{
  //
  // default contstructor: do nothing
  //
}

//____________________________________________________________________________________//

AliAnalysisTaskTpcSpectra::AliAnalysisTaskTpcSpectra(const char* name) : AliAnalysisTaskSE(name), fAOD(0), fOutputList(0), fEventCut(0),
                                                                         //fAODTrackCuts(0),
                                                                         //fAODTrackCutsWoDCA(0),
                                                                         fPIDResponse(0),
                                                                         fMultSel(0),
                                                                         fHistZv(0),
                                                                         fHistdEdxData(0),
                                                                         fHistTof(0),
                                                                         fHistCentBefEvSel(0),
                                                                         fHistCent(0),
                                                                         fHistoMult(0),
                                                                         fHistMuonGen(0),
                                                                         fHistMuonReco(0),
                                                                         fHistMuonRecoTOF(0),
                                                                         fHistPionGen(0),
                                                                         fHistPionReco(0),
                                                                         fHistPionRecoTOF(0),
                                                                         fHistKaonGen(0),
                                                                         fHistKaonReco(0),
                                                                         fHistKaonRecoTOF(0),
                                                                         fHistProtGen(0),
                                                                         fHistProtReco(0),
                                                                         fHistProtRecoTOF(0),
                                                                         fHistAntiMuonGen(0),
                                                                         fHistAntiMuonReco(0),
                                                                         fHistAntiMuonRecoTOF(0),
                                                                         fHistAntiPionGen(0),
                                                                         fHistAntiPionReco(0),
                                                                         fHistAntiPionRecoTOF(0),
                                                                         fHistAntiKaonGen(0),
                                                                         fHistAntiKaonReco(0),
                                                                         fHistAntiKaonRecoTOF(0),
                                                                         fHistAntiProtGen(0),
                                                                         fHistAntiProtReco(0),
                                                                         fHistAntiProtRecoTOF(0),
                                                                         fMuonDeDxCent(0),
                                                                         fPionDeDxCent(0),
                                                                         fKaonDeDxCent(0),
                                                                         fProtonDeDxCent(0),
                                                                         fAntiMuonDeDxCent(0),
                                                                         fAntiPionDeDxCent(0),
                                                                         fAntiKaonDeDxCent(0),
                                                                         fAntiProtonDeDxCent(0),
                                                                         fDCAPion(0),
                                                                         fDCAAntiPion(0),
                                                                         fDCAKaon(0),
                                                                         fDCAAntiKaon(0),
                                                                         fDCAProton(0),
                                                                         fDCAAntiProton(0),
                                                                         fDCAPionMC(0),
                                                                         fDCAAntiPionMC(0),
                                                                         fDCAKaonMC(0),
                                                                         fDCAAntiKaonMC(0),
                                                                         fDCAProtonMC(0),
                                                                         fDCAAntiProtonMC(0)
{
  //
  // main constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________________________//

AliAnalysisTaskTpcSpectra::~AliAnalysisTaskTpcSpectra()
{
  //
  // destructor
  //
  if (fOutputList) {
    delete fOutputList;
  }
}

//____________________________________________________________________________________//

void AliAnalysisTaskTpcSpectra::UserCreateOutputObjects()
{
  //
  // create the output objects
  //
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  fEventCut.AddQAplotsToList(fOutputList);
  //
  // pt-bins
  //
  const Int_t nPtBins = 59;
  Double_t binsPt[nPtBins + 1] = {0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35,
                                  0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
                                  0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
                                  1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40,
                                  3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00,
                                  9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 18.00, 20.00};
  //
  // create QA histograms
  //
  fHistZv = new TH1F("fHistZv", "fHistZv; vertex z (cm)", 200, -40, 40);
  fHistdEdxData = new TH2F("fHistdEdxData", "fHistdEdxData; p/z (GeV/c); TPC signal", 500, -10, 10, 250, 0, 500);
  fHistTof = new TH2F("fHistTof", "all paritlces TOF; P(Gev/c); beta", 1000, -10, 10, 1000, 0, 1.5); // all particles TOF quality assurance
  fHistCentBefEvSel = new TH1F("fHistCentBefEvSel", "fHistCentBefEvSel", 110, -0.01, 0.1);
  fHistCent = new TH1F("fHistCent", "fHistCent", 110, -0.01, 0.1);
  fHistoMult = new TH1F("fHistoMult", "fHistoMult;  z", 200, 0, 200);
  //
  // create histograms for efficiency calculation
  //
  fHistMuonGen = new TH1F("fHistMuonGen", "Generated Muon #phi vs p_{t}", nPtBins, binsPt);
  fHistMuonReco = new TH1F("fHistMuonReco", "Reconstructed Muon #phi vs p_{t}", nPtBins, binsPt);
  fHistMuonRecoTOF = new TH1F("fHistMuonRecoTOF", "Reconstructed Muon #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistPionGen = new TH1F("fHistPionGen", "Generated Pion #phi vs p_{t}", nPtBins, binsPt);
  fHistPionReco = new TH1F("fHistPionReco", "Reconstructed Pion #phi vs p_{t}", nPtBins, binsPt);
  fHistPionRecoTOF = new TH1F("fHistPionRecoTOF", "Reconstructed Pion #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistKaonGen = new TH1F("fHistKaonGen", "Generated Kaon #phi vs p_{t}", nPtBins, binsPt);
  fHistKaonReco = new TH1F("fHistKaonReco", "Reconstructed Kaon #phi vs p_{t}", nPtBins, binsPt);
  fHistKaonRecoTOF = new TH1F("fHistKaonRecoTOF", "Reconstructed Kaon #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistProtGen = new TH1F("fHistProtGen", "Generated Proton #phi vs p_{t}", nPtBins, binsPt);
  fHistProtReco = new TH1F("fHistProtReco", "Reconstructed Proton #phi vs p_{t}", nPtBins, binsPt);
  fHistProtRecoTOF = new TH1F("fHistProtRecoTOF", "Reconstructed Proton #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistAntiMuonGen = new TH1F("fHistAntiMuonGen", "Generated AntiMuon #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiMuonReco = new TH1F("fHistAntiMuonReco", "Reconstructed AntiMuon #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiMuonRecoTOF = new TH1F("fHistAntiMuonRecoTOF", "Reconstructed AntiMuon #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistAntiPionGen = new TH1F("fHistAntiPionGen", "Generated AntiPion #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiPionReco = new TH1F("fHistAntiPionReco", "Reconstructed AntiPion #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiPionRecoTOF = new TH1F("fHistAntiPionRecoTOF", "Reconstructed AntiPion #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistAntiKaonGen = new TH1F("fHistAntiKaonGen", "Generated AntiKaon #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiKaonReco = new TH1F("fHistAntiKaonReco", "Reconstructed AntiKaon #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiKaonRecoTOF = new TH1F("fHistAntiKaonRecoTOF", "Reconstructed AntiKaon #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistAntiProtGen = new TH1F("fHistAntiProtGen", "Generated AntiProton #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiProtReco = new TH1F("fHistAntiProtReco", "Reconstructed AntiProton #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiProtRecoTOF = new TH1F("fHistAntiProtRecoTOF", "Reconstructed AntiProton #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  // create raw spectra histograms
  //
  fMuonDeDxCent = new TH3F("fMuonDeDxCent", "fMuonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -0.05, 0.15);
  fPionDeDxCent = new TH3F("fPionDeDxCent", "fPionDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -0.05, 0.15);
  fKaonDeDxCent = new TH3F("fKaonDeDxCent", "fKaonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -0.05, 0.15);
  fProtonDeDxCent = new TH3F("fProtonDeDxCent", "fProtonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -0.05, 0.15);
  fMuonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fPionDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fKaonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fProtonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  //
  fAntiMuonDeDxCent = new TH3F("fAntiMuonDeDxCent", "fAntiMuonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -0.05, 0.15);
  fAntiPionDeDxCent = new TH3F("fAntiPionDeDxCent", "fAntiPionDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -0.05, 0.15);
  fAntiKaonDeDxCent = new TH3F("fAntiKaonDeDxCent", "fAntiKaonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -0.05, 0.15);
  fAntiProtonDeDxCent = new TH3F("fAntiProtonDeDxCent", "fAntiProtonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -0.05, 0.15);
  fAntiMuonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fAntiPionDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fAntiKaonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fAntiProtonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  //
  // create feed-down correction histograms
  //
  fDCAPion = new TH3F("fDCAPion", "fDCAPion", nPtBins, 0, 20, 600, -3, 3, 200, -0.05, 0.15);
  fDCAAntiPion = new TH3F("fDCAAntiPion", "fDCAAntiPion", nPtBins, 0, 20, 600, -3, 3, 200, -0.05, 0.15);
  fDCAKaon = new TH3F("fDCAKaon", "fDCAKaon", nPtBins, 0, 20, 600, -3, 3, 200, -0.05, 0.15);
  fDCAAntiKaon = new TH3F("fDCAAntiKaon", "fDCAAntiKaon", nPtBins, 0, 20, 600, -3, 3, 200, -0.05, 0.15);
  fDCAProton = new TH3F("fDCAProton", "fDCAProton", nPtBins, 0, 20, 600, -3, 3, 200, -0.05, 0.15);
  fDCAAntiProton = new TH3F("fDCAAntiProton", "fDCAAntiProton", nPtBins, 0, 20, 600, -3, 3, 200, -0.05, 0.150);
  fDCAPion->GetXaxis()->Set(nPtBins, binsPt);
  fDCAAntiPion->GetXaxis()->Set(nPtBins, binsPt);
  fDCAKaon->GetXaxis()->Set(nPtBins, binsPt);
  fDCAAntiKaon->GetXaxis()->Set(nPtBins, binsPt);
  fDCAProton->GetXaxis()->Set(nPtBins, binsPt);
  fDCAAntiProton->GetXaxis()->Set(nPtBins, binsPt);
  //
  // create feed-down correction histograms for MC
  //
  fDCAPionMC = new TH3F("fDCAPionMC", "fDCAPionMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5);
  fDCAAntiPionMC = new TH3F("fDCAAntiPionMC", "fDCAAntiPionMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5);
  fDCAKaonMC = new TH3F("fDCAKaonMC", "fDCAKaonMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5);
  fDCAAntiKaonMC = new TH3F("fDCAAntiKaonMC", "fDCAAntiKaonMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5);
  fDCAProtonMC = new TH3F("fDCAProtonMC", "fDCAProtonMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5);
  fDCAAntiProtonMC = new TH3F("fDCAAntiProtonMC", "fDCAAntiProtonMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5);
  fDCAPionMC->GetXaxis()->Set(nPtBins, binsPt);
  fDCAAntiPionMC->GetXaxis()->Set(nPtBins, binsPt);
  fDCAKaonMC->GetXaxis()->Set(nPtBins, binsPt);
  fDCAAntiKaonMC->GetXaxis()->Set(nPtBins, binsPt);
  fDCAProtonMC->GetXaxis()->Set(nPtBins, binsPt);
  fDCAAntiProtonMC->GetXaxis()->Set(nPtBins, binsPt);
  //add histograms to output
  fOutputList->Add(fHistZv);
  fOutputList->Add(fHistdEdxData);
  fOutputList->Add(fHistTof);
  fOutputList->Add(fHistCentBefEvSel);
  fOutputList->Add(fHistCent);
  fOutputList->Add(fHistoMult);
  //
  fOutputList->Add(fHistMuonGen);
  fOutputList->Add(fHistMuonReco);
  fOutputList->Add(fHistMuonRecoTOF);
  //
  fOutputList->Add(fHistPionGen);
  fOutputList->Add(fHistPionReco);
  fOutputList->Add(fHistPionRecoTOF);
  //
  fOutputList->Add(fHistKaonGen);
  fOutputList->Add(fHistKaonReco);
  fOutputList->Add(fHistKaonRecoTOF);
  //
  fOutputList->Add(fHistProtGen);
  fOutputList->Add(fHistProtReco);
  fOutputList->Add(fHistProtRecoTOF);
  //
  fOutputList->Add(fHistAntiMuonGen);
  fOutputList->Add(fHistAntiMuonReco);
  fOutputList->Add(fHistAntiMuonRecoTOF);
  //
  fOutputList->Add(fHistAntiPionGen);
  fOutputList->Add(fHistAntiPionReco);
  fOutputList->Add(fHistAntiPionRecoTOF);
  //
  fOutputList->Add(fHistAntiKaonGen);
  fOutputList->Add(fHistAntiKaonReco);
  fOutputList->Add(fHistAntiKaonRecoTOF);
  //
  fOutputList->Add(fHistAntiProtGen);
  fOutputList->Add(fHistAntiProtReco);
  fOutputList->Add(fHistAntiProtRecoTOF);
  //
  fOutputList->Add(fMuonDeDxCent);
  fOutputList->Add(fPionDeDxCent);
  fOutputList->Add(fKaonDeDxCent);
  fOutputList->Add(fProtonDeDxCent);
  //
  fOutputList->Add(fAntiMuonDeDxCent);
  fOutputList->Add(fAntiPionDeDxCent);
  fOutputList->Add(fAntiKaonDeDxCent);
  fOutputList->Add(fAntiProtonDeDxCent);
  //
  fOutputList->Add(fDCAPion);
  fOutputList->Add(fDCAAntiPion);
  fOutputList->Add(fDCAKaon);
  fOutputList->Add(fDCAAntiKaon);
  fOutputList->Add(fDCAProton);
  fOutputList->Add(fDCAAntiProton);
  //
  fOutputList->Add(fDCAPionMC);
  fOutputList->Add(fDCAAntiPionMC);
  fOutputList->Add(fDCAKaonMC);
  fOutputList->Add(fDCAAntiKaonMC);
  fOutputList->Add(fDCAProtonMC);
  fOutputList->Add(fDCAAntiProtonMC);

  //
  PostData(1, fOutputList);
}

//____________________________________________________________________________________//
void AliAnalysisTaskTpcSpectra::UserExec(Option_t*)
{
  //
  // loop over events
  //
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();  
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    if (inputHandler)
      fPIDResponse = inputHandler->GetPIDResponse();
  }
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD)
    return;
  bool EvtAcc = fEventCut.AcceptEvent(fAOD);
  if (!fEventCut.PassedCut(AliEventCuts::kTrigger))
    return;
  //
  //track multiplicity at mid-rapidity (counts number of tracks we reconstruct)
  //
  AliAODHeader* aodHeader = (AliAODHeader*)fAOD->GetHeader();
  fHistoMult->Fill(aodHeader->GetRefMultiplicity());
  //
  // centrality
  //
  Float_t centrality = -999;
  fMultSel = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if (!fMultSel) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  } else {
    AliDebug(2, "Estimating centrality");
    centrality = fMultSel->GetMultiplicityPercentile("V0M", kFALSE); //Event selection is embedded in the Multiplicity estimator so that the Multiplicit                                                                      y percentiles are well defined and refer to the same sample
  }
  //Fill centrality histogram before event selection
  fHistCentBefEvSel->Fill(centrality);
  if (!EvtAcc)
    return;
  //
  // primary vertex
  //
  //const AliVVertex *vertex = fAOD->GetPrimaryVertexTracks();
  AliAODVertex* vertex = (AliAODVertex*)fAOD->GetPrimaryVertexTracks();
  if (!vertex || (vertex && vertex->GetNContributors() < 1)) {


    // SPD vertex

    vertex = fAOD->GetPrimaryVertexSPD();
    if (!vertex)
      return;
    if (vertex->GetNContributors() < 1)
      vertex = 0x0;
  }
  if (!vertex)
    return;
  fHistZv->Fill(vertex->GetZ());
  if (TMath::Abs(vertex->GetZ()) > 10.0)
    return; // remove events with a vertex which is more than 10cm away
  //Fill centrality histogram after event selection
  fHistCent->Fill(centrality);
  //
  // MONTE CARLO -- GENERATED PARTICLES
  //
  //Bool_t mcTrue = kFALSE;
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  //if (eventHandler) mcTrue = kTRUE;
  //
  AliMCEvent* mcEvent = 0x0;
  if (eventHandler)
    mcEvent = MCEvent();
    else Printf("Didnot find event handler");

  if (mcEvent) {
    //
    // loop over generated particles
    //
    for (Int_t i = 0; i < mcEvent->GetNumberOfTracks(); i++) { // start loop generated particles
//AliAODMCParticle* AODParticle = (AliAODMCParticle*)mcEvent->GetTrack(abs(track->GetLabel()));
    AliAODMCParticle* trackMC = (AliAODMCParticle*)mcEvent-> GetTrack(i);  
    //  TParticle* trackMC = ((AliMCParticle*)mcEvent->GetTrack(i))->Particle();
      if (!trackMC)
        continue;
      Int_t mcCode = trackMC->GetPdgCode();

      if (trackMC->IsPhysicalPrimary() && TMath::Abs(trackMC->Y()) < 0.5) {
        switch (mcCode) {
          case 13:
            fHistMuonGen->Fill(trackMC->Pt());
            break; //muon
          case 211: fHistPionGen->Fill(trackMC->Pt()); break; //pion
          case 321: fHistKaonGen->Fill(trackMC->Pt()); break; //kaon
          case 2212: fHistProtGen->Fill(trackMC->Pt()); break; //proton
          case -13: fHistAntiMuonGen->Fill(trackMC->Pt()); break; //Antimuon
          case -211: fHistAntiPionGen->Fill(trackMC->Pt()); break; //Antipion
          case -321: fHistAntiKaonGen->Fill(trackMC->Pt()); break; //Antikaon
          case -2212: fHistAntiProtGen->Fill(trackMC->Pt()); break; //antiproton
        }
      }
    } // end loop generated particles
  } else
    Printf("*****Didnot find MC events*****");
  //
  // RECONSTRUCTED TRACKS
  //
  Int_t jTracks = fAOD->GetNumberOfTracks();
  Float_t dcaxy, dcaz;
  Float_t Prod = 0;
  for (Int_t j = 0; j < jTracks; j++) { // start loop over tracks
    //
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
    if (!track) {
      continue;
    }
    if (TMath::Abs(track->Eta()) > 0.8) {
      continue;
    }
    if (!track->TestFilterBit(32)) {
      continue;
    }
    //if (!fAODTrackCutsWoDCA->AcceptTrack(track)) continue; // check if track passes the cuts
    //Get the DCAxy and DCAz
    track->GetImpactParameters(dcaxy, dcaz);
    //
    // calculate rapidities and kinematics
    //
    Double_t pvec[3];
    track->GetPxPyPz(pvec);
    Double_t energyMuon = TMath::Sqrt(track->GetP() * track->GetP() + AliPID::ParticleMass(AliPID::kMuon) * AliPID::ParticleMass(AliPID::kMuon));   
    Double_t energyPion = TMath::Sqrt(track->GetP() * track->GetP() + AliPID::ParticleMass(AliPID::kPion) * AliPID::ParticleMass(AliPID::kPion));
    Double_t energyKaon = TMath::Sqrt(track->GetP() * track->GetP() + AliPID::ParticleMass(AliPID::kKaon) * AliPID::ParticleMass(AliPID::kKaon));
    Double_t energyProton = TMath::Sqrt(track->GetP() * track->GetP() + AliPID::ParticleMass(AliPID::kProton) * AliPID::ParticleMass(AliPID::kProton));
    //
      const Double_t abspl = TMath::Abs(pvec[2]); 
    Double_t rapMuon = energyMuon > abspl ? 0.5 * TMath::Log((energyMuon + pvec[2]) / (energyMuon - pvec[2])) : 999.;
      Printf("rapMuon%f",rapMuon);
    Double_t rapPion = energyPion > abspl ? 0.5 * TMath::Log((energyPion + pvec[2]) / (energyPion - pvec[2])) : 999.;
    Double_t rapKaon = energyKaon > abspl ? 0.5 * TMath::Log((energyKaon + pvec[2]) / (energyKaon - pvec[2])) : 999.;
    Double_t rapProton = energyProton > abspl ? 0.5 * TMath::Log((energyProton + pvec[2]) / (energyProton - pvec[2])) : 999.;
    if (track->Charge() > 0) {
      if (TMath::Abs(rapPion) < 0.5 && GetCombinedSigma(track, AliPID::kPion) < 2)
        fDCAPion->Fill(track->Pt(), dcaxy, centrality);
      if (TMath::Abs(rapKaon) < 0.5 && GetCombinedSigma(track, AliPID::kKaon) < 2)
        fDCAKaon->Fill(track->Pt(), dcaxy, centrality);
      if (TMath::Abs(rapProton) < 0.5 && GetCombinedSigma(track, AliPID::kProton) < 2)
        fDCAProton->Fill(track->Pt(), dcaxy, centrality);
    } 
    else if (track->Charge() < 0) {
      if (TMath::Abs(rapPion) < 0.5 && GetCombinedSigma(track, AliPID::kPion) < 2)
        fDCAAntiPion->Fill(track->Pt(), dcaxy, centrality);
      if (TMath::Abs(rapKaon) < 0.5 && GetCombinedSigma(track, AliPID::kKaon) < 2)
        fDCAAntiKaon->Fill(track->Pt(), dcaxy, centrality);
      if (TMath::Abs(rapProton) < 0.5 && GetCombinedSigma(track, AliPID::kProton) < 2)
        fDCAAntiProton->Fill(track->Pt(), dcaxy, centrality);
    }
    if (mcEvent) {
      AliAODMCParticle* AODParticle = (AliAODMCParticle*)mcEvent->GetTrack(abs(track->GetLabel()));
      //TParticle* trackReco = ((AliMCParticle*)mcEvent->GetTrack(abs(track->GetLabel())))->Particle();

      Int_t recoCode = AODParticle->GetPdgCode();
      if (TMath::Abs(AODParticle->Y()) < 0.5) {
        Prod = 3;
        if (AODParticle->IsPhysicalPrimary())
          Prod = 0;
        else if (AODParticle->IsSecondaryFromWeakDecay())
          Prod = 1;
        else if (AODParticle->IsSecondaryFromMaterial())
          Prod = 2;
        //
        if (TMath::Abs(recoCode) == AliPID::ParticleCode(AliPID::kPion)) {
          if (recoCode > 0)
            fDCAPionMC->Fill(track->Pt(), dcaxy, Prod);
          else
            fDCAAntiPionMC->Fill(track->Pt(), dcaxy, Prod);
        } else if (TMath::Abs(recoCode) == AliPID::ParticleCode(AliPID::kKaon)) {
          if (recoCode > 0)
            fDCAKaonMC->Fill(track->Pt(), dcaxy, Prod);
          else
            fDCAAntiKaonMC->Fill(track->Pt(), dcaxy, Prod);
        } else if (TMath::Abs(recoCode) == AliPID::ParticleCode(AliPID::kProton)) {
          if (recoCode > 0)
            fDCAProtonMC->Fill(track->Pt(), dcaxy, Prod);
          else
            fDCAAntiProtonMC->Fill(track->Pt(), dcaxy, Prod);
        }
      }
    } else
    Printf("*****Didnot find MC events*****");
    //if (!fAODTrackCuts->AcceptTrack(track)) continue; // check if track passes the cuts
   
    Printf(" %f\n", track->P());
    //if (!track->GetInnerParam()) continue;            // check if track is a proper TPC track
    //
    Double_t ptot = track->P();       // momentum for dEdx determination
    Double_t sign = track->GetSign(); // charge
    Printf("%f\n", sign);
    //
    // fill the QA dE/Dx histograms
    //
    //Printf("************************TPC SIGNAL***********************");
    //Printf(" %f",track->GetTPCsignal());
    if (track->GetTPCsignalN() > 80) {
      Printf("Filling %f %f\n", ptot * sign, track->GetTPCsignal());
      fHistdEdxData->Fill(ptot * sign, track->GetTPCsignal());
    }
    //
    // fill the raw spectra
    //
    if (track->Charge() > 0) {
      if (TMath::Abs(rapMuon) < 0.5)
        fMuonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon), centrality);
      if (TMath::Abs(rapPion) < 0.5)
        fPionDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion), centrality);
      if (TMath::Abs(rapKaon) < 0.5)
        fKaonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon), centrality);
      if (TMath::Abs(rapProton) < 0.5)
        fProtonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton), centrality);
    } else if (track->Charge() < 0) {
      if (TMath::Abs(rapMuon) < 0.5)
        fAntiMuonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon), centrality);
      if (TMath::Abs(rapPion) < 0.5)
        fAntiPionDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion), centrality);
      if (TMath::Abs(rapKaon) < 0.5)
        fAntiKaonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon), centrality);
      if (TMath::Abs(rapProton) < 0.5)
        fAntiProtonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton), centrality);
    }
    //
    //Check TOF status
    //
    UInt_t status = 0;
    status = track->GetStatus();
    //
    Bool_t hasTOFout = status & AliAODTrack::kTOFout;
    Bool_t hasTOFtime = status & AliAODTrack::kTIME;
    Bool_t hasTOF = kFALSE;
    if (hasTOFout && hasTOFtime)
      hasTOF = kTRUE;
    Float_t length = track->GetIntegratedLength();
    if (length < 350.)
      hasTOF = kFALSE;
    //
    // access mcTruth for efficiency calculation
    //
    if (track->GetTPCsignal() < 10.0) continue; // make sure that PID efficiency is taken into account
    if (mcEvent) {
      AliAODMCParticle* AODParticle = (AliAODMCParticle*)mcEvent->GetTrack(abs(track->GetLabel()));
      //TParticle* trackReco = ((AliMCParticle*)mcEvent->GetTrack(abs(track->GetLabel())))->Particle();
      Int_t recoCode = AODParticle->GetPdgCode();
      //Without TOF
      //if (mcEvent->IsPhysicalPrimary(abs(track->GetLabel())) && abs(trackReco->Y()) < 0.5)
      if ((AODParticle->IsPhysicalPrimary()) && abs(AODParticle->Y()) < 0.5) {
        switch (recoCode) {
          case 13:
            fHistMuonReco->Fill(AODParticle->Pt());
            break; //muon
          case 211:
            fHistPionReco->Fill(AODParticle->Pt());
            break; //pion
          case 321:
            fHistKaonReco->Fill(AODParticle->Pt());
            break; //kaon
          case 2212:
            fHistProtReco->Fill(AODParticle->Pt());
            break; //proton
          case -13:
            fHistAntiMuonReco->Fill(AODParticle->Pt());
            break; //Antimuon
          case -211:
            fHistAntiPionReco->Fill(AODParticle->Pt());
            break; //Antipion
          case -321:
            fHistAntiKaonReco->Fill(AODParticle->Pt());
            break; //Antikaon
          case -2212:
            fHistAntiProtReco->Fill(AODParticle->Pt());
            break; //Antiproton
        }
      }
      //With TOF
      if (AODParticle->IsPhysicalPrimary() && abs(AODParticle->Y()) < 0.5 && track->GetTPCncls() > 70 && hasTOF) {
        switch (recoCode) {
          case 13:
            fHistMuonRecoTOF->Fill(AODParticle->Pt());
            break; //muon
          case 211:
            fHistPionRecoTOF->Fill(AODParticle->Pt());
            break; //pion
          case 321:
            fHistKaonRecoTOF->Fill(AODParticle->Pt());
            break; //kaon
          case 2212:
            fHistProtRecoTOF->Fill(AODParticle->Pt());
            break; //proton
          case -13:
            fHistAntiMuonRecoTOF->Fill(AODParticle->Pt());
            break; //Antimuon
          case -211:
            fHistAntiPionRecoTOF->Fill(AODParticle->Pt());
            break; //Antipion
          case -321:
            fHistAntiKaonRecoTOF->Fill(AODParticle->Pt());
            break; //Antikaon
          case -2212:
            fHistAntiProtRecoTOF->Fill(AODParticle->Pt());
            break; //Antiproton
        }
      }
    } // end access to MC truth

  } // end track loop
  //
  // post the data and end the event loop
  //
  PostData(1, fOutputList);
}

//____________________________________________________________________________________//

void AliAnalysisTaskTpcSpectra::Terminate(Option_t*)
{
}

//____________________________________________________________________________________//
