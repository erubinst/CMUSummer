// -*- C++ -*-
#if !defined(jetHists_H)
#define jetHists_H

#include <iostream>
#include <TH1F.h>
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "nTupleAnalysis/baseClasses/interface/jetData.h"
#include "nTupleAnalysis/baseClasses/interface/fourVectorHists.h"

namespace nTupleAnalysis {

  class jetHists {
  public:
    TFileDirectory dir;
    
    fourVectorHists* v;
    TH1F* cleanmask;
    TH1F* deepB;
    TH1F* CSVv2;
    TH1F* deepFlavB;

    jetHists(std::string, fwlite::TFileService&, std::string title = "");
    void Fill(std::shared_ptr<jet>&, float);
    ~jetHists(); 

  };

}
#endif // jetHists_H
