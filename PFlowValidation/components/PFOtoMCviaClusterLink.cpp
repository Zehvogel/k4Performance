/*
 * PFOtoMCviaClusterLink: PFOs in signature, Reco↔MC links used to match,
 * histograms saved via THistSvc (no TEfficiency, no TVector3).
 */

#include "k4FWCore/Consumer.h"
#include "k4FWCore/DataHandle.h"

#include "Gaudi/Property.h"
#include "GaudiKernel/ITHistSvc.h"

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/RecoMCParticleLinkCollection.h"

#include "podio/ObjectID.h"

#include "TH1F.h"
#include "TH2F.h"

#include <map>
#include <vector>
#include <string>
#include <cmath>

namespace {
  // Works for edm4hep::Vector3f and Vector3d
  template <typename V3>
  inline float safeCosTheta(const V3& p) {
    const double px = static_cast<double>(p.x);
    const double py = static_cast<double>(p.y);
    const double pz = static_cast<double>(p.z);
    const double mag = std::sqrt(px*px + py*py + pz*pz);
    return mag > 0. ? float(pz / mag) : 1.f;
  }
  struct ObjectIDLess {
    bool operator()(const podio::ObjectID& a, const podio::ObjectID& b) const noexcept {
      if (a.collectionID != b.collectionID) return a.collectionID < b.collectionID;
      return a.index < b.index;
    }
  };
  inline int pidBin(int pdg) {
    switch (pdg) { case 22: return 0; case 11: return 1; case -11: return 2; case 211: return 3; case -211: return 4; default: return 5; }
  }
}

struct PFOtoMCviaClusterLink final
  : k4FWCore::Consumer<void(const edm4hep::ReconstructedParticleCollection&,
                            const edm4hep::MCParticleCollection&,
                            const edm4hep::ClusterCollection&,
                            const edm4hep::RecoMCParticleLinkCollection&)> {
  PFOtoMCviaClusterLink(const std::string& name, ISvcLocator* svcLoc)
    : Consumer(name, svcLoc, {
        KeyValues("InputPFOs", {"PandoraPFOs"}),
        KeyValues("InputMCParticles", {"MCParticles"}),
        KeyValues("InputClusters", {"PandoraClusters"}),
        KeyValues("InputRecoMC", {"MCTruthRecoLink"})
      }) {}

  // Properties
  Gaudi::Property<std::string> m_histPath{this,"HistPath","/MC/Cluster",
      "THistSvc directory (must start with the stream name, e.g. '/MC/...')."};
  Gaudi::Property<int>    m_binsE{this,"BinsE",64,"Bins in energy [GeV]."};
  Gaudi::Property<double> m_eMin{this,"EMin",0.,"Min E [GeV]."};
  Gaudi::Property<double> m_eMax{this,"EMax",128.,"Max E [GeV]."};
  Gaudi::Property<int>    m_binsCos{this,"BinsCosT",30,"Bins in cos(theta)."};

  // Services & hists
  mutable SmartIF<ITHistSvc> m_histSvc;
  mutable bool m_booked=false;

  // Simple spectra
  mutable TH1F *h_mcE=nullptr, *h_pfoE=nullptr, *h_mcCos=nullptr, *h_pfoCos=nullptr;
  // Response (Reco vs MC)
  mutable TH2F *h_respE=nullptr, *h_respCos=nullptr;
  // Efficiency as num/den
  mutable TH1F *h_effE_num=nullptr, *h_effE_den=nullptr;
  mutable TH1F *h_effCos_num=nullptr, *h_effCos_den=nullptr;

  void bookIfNeeded() const {
    if (m_booked) return;
    m_histSvc = service("THistSvc");
    if (!m_histSvc) { error()<<"THistSvc not found"<<endmsg; return; }
    auto reg = [&](TH1* h){ return m_histSvc->regHist(m_histPath.value()+"/"+h->GetName(), h).isSuccess(); };

    h_mcE   = new TH1F("mcE","MC E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_pfoE  = new TH1F("pfoE","PFO E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_mcCos = new TH1F("mcCos","MC cos#theta;cos#theta;Entries", m_binsCos, -1, 1);
    h_pfoCos= new TH1F("pfoCos","PFO cos#theta;cos#theta;Entries", m_binsCos, -1, 1);

    h_respE   = new TH2F("respE","Reco E vs MC E;MC E [GeV];Reco E [GeV]", m_binsE, m_eMin, m_eMax, m_binsE, m_eMin, m_eMax);
    h_respCos = new TH2F("respCos","Reco cos vs MC cos;MC cos#theta;Reco cos#theta", m_binsCos, -1, 1, m_binsCos, -1, 1);

    h_effE_num   = new TH1F("effE_num","Efficiency numerator vs E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_effE_den   = new TH1F("effE_den","Efficiency denominator vs E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_effCos_num = new TH1F("effCos_num","Efficiency numerator vs cos#theta;cos#theta;Entries", m_binsCos, -1, 1);
    h_effCos_den = new TH1F("effCos_den","Efficiency denominator vs cos#theta;cos#theta;Entries", m_binsCos, -1, 1);

    bool ok=true;
    ok&=reg(h_mcE); ok&=reg(h_pfoE); ok&=reg(h_mcCos); ok&=reg(h_pfoCos);
    ok&=reg(h_respE); ok&=reg(h_respCos);
    ok&=reg(h_effE_num); ok&=reg(h_effE_den); ok&=reg(h_effCos_num); ok&=reg(h_effCos_den);
    if (!ok) error()<<"Failed to register some histograms (check HistPath stream)."<<endmsg;

    m_booked=true;
  }

  void operator()(const edm4hep::ReconstructedParticleCollection& pfos, const edm4hep::MCParticleCollection& mcps, const edm4hep::ClusterCollection& clus, const edm4hep::RecoMCParticleLinkCollection& links) const override {
    bookIfNeeded(); if (!m_booked) return;

    // MC maps
    std::map<podio::ObjectID, size_t, ObjectIDLess> mcIndex;
    for (size_t i=0;i<mcps.size();++i) mcIndex[mcps[i].getObjectID()] = i;

    // Fill MC denoms/spectra once per event
    for (const auto& mc : mcps) {
      h_mcE->Fill(mc.getEnergy());
      h_mcCos->Fill(safeCosTheta(mc.getMomentum()));
      h_effE_den->Fill(mc.getEnergy());
      h_effCos_den->Fill(safeCosTheta(mc.getMomentum()));
    }

    // Group RecoMC links by reco id
    std::map<podio::ObjectID, std::vector<const edm4hep::RecoMCParticleLink*>, ObjectIDLess> byReco;
    for (const auto& L : links) {
      if (!L.getFrom().isAvailable() || !L.getTo().isAvailable()) continue;  // from: Reco, to: MC
      byReco[L.getFrom().getObjectID()].push_back(&L);
    }

    // Loop over PFOs
    for (const auto& rp : pfos) {
      const float eRec = rp.getEnergy();
      const float cRec = safeCosTheta(rp.getMomentum());
      h_pfoE->Fill(eRec);
      h_pfoCos->Fill(cRec);

      // pick best linked MC by weight
      int bestIdx=-1; float bestW=-1.f;
      auto itL = byReco.find(rp.getObjectID());
      if (itL != byReco.end()) {
        for (auto Lptr : itL->second) {
          const auto mid = Lptr->getTo().getObjectID();
          auto itMC = mcIndex.find(mid);
          if (itMC==mcIndex.end()) continue;
          const float w = Lptr->getWeight();
          if (w > bestW) { bestW = w; bestIdx = (int)itMC->second; }
        }
      }

      if (bestIdx>=0) {
        const auto& mc = mcps[bestIdx];
        const float eMC = mc.getEnergy();
        const float cMC = safeCosTheta(mc.getMomentum());
        h_respE->Fill(eMC, eRec);
        h_respCos->Fill(cMC, cRec);
        h_effE_num->Fill(eMC);
        h_effCos_num->Fill(cMC);
      }
    }
  }
};

DECLARE_COMPONENT(PFOtoMCviaClusterLink)
