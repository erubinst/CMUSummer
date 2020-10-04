#include "AnaExample/nTupleAnalysis/interface/eventHists.h"

using namespace nTupleAnalysis;

eventHists::eventHists(std::string name, fwlite::TFileService& fs, bool isMC) {
  std::cout << "Initialize >> eventHists: " << name << std::endl;
  dir = fs.mkdir(name);

  //
  // Object Level
  //
  nAllJets = dir.make<TH1F>("nAllJets", (name+"/nAllJets; Number of Jets (no selection); Entries").c_str(),  16,-0.5,15.5);
  nSelJets = dir.make<TH1F>("nSelJets", (name+"/nSelJets; Number of Selected Jets; Entries").c_str(),  16,-0.5,15.5);
  allJets = new jetHists(name+"/allJets", fs, "All Jets");
  selJets = new jetHists(name+"/selJets", fs, "Selected Jets");
    
  nAllMuons = dir.make<TH1F>("nAllMuons", (name+"/nAllMuons; Number of Muons (no selection); Entries").c_str(),  6,-0.5,5.5);
  nIsoMuons = dir.make<TH1F>("nIsoMuons", (name+"/nIsoMuons; Number of Prompt Muons; Entries").c_str(),  6,-0.5,5.5);
  allMuons = new muonHists(name+"/allMuons", fs, "All Muons");
  isoMuons = new muonHists(name+"/isoMuons", fs, "Prompt Muons");

} 

void eventHists::Fill(eventData* event){
  //
  // Object Level
  //
  nAllJets->Fill(event->allJets.size(), event->weight);
  nSelJets->Fill(event->selJets.size(), event->weight);
  for(auto &jet: event->allJets) allJets->Fill(jet, event->weight);
  for(auto &jet: event->selJets) selJets->Fill(jet, event->weight);

  nAllMuons->Fill(event->allMuons.size(), event->weight);
  nIsoMuons->Fill(event->isoMuons.size(), event->weight);
  for(auto &muon: event->allMuons) allMuons->Fill(muon, event->weight);
  for(auto &muon: event->isoMuons) isoMuons->Fill(muon, event->weight);

  return;
}


eventHists::~eventHists(){} 

