/*
 * PFOtoMCviaClusterLink: PFOs in signature, Reco↔MC links used to match,
 * histograms saved via THistSvc (no TEfficiency, no TVector3).
 *
 * Adds split photon efficiencies (converted/unconverted) vs E and vs cosθ.
 */

#include "k4FWCore/Consumer.h"

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
#include <algorithm>

namespace {
  // Works for edm4hep::Vector3f and Vector3d (x,y,z public members)
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
  inline int pidCategory(int pdg) {
    if (pdg==13 || pdg==-13) return 0;      // mu
    if (pdg==22)             return 1;      // gamma
    if (pdg==11 || pdg==-11) return 2;      // e
    if (pdg==130 || pdg==2112) return 3;    // neutral hadron (K0L,n)
    if (pdg==211 || pdg==-211 || pdg==2212 || pdg==-2212) return 4; // charged hadron (pi, p)
    return 5;
  }

  inline bool isElectronPDG(int pdg){ return pdg==11 || pdg==-11; }

  // Truth definition of 'converted photon':
  // photon with >=2 e+e- daughters AND vertex R within [rMinConv, rMaxConv] (configurable)
  template <typename V3>
  inline bool isConvertedPhoton(const edm4hep::MCParticle& g, float rMinConv, float rMaxConv) {
    if (g.getPDG() != 22) return false;
    int nEle = 0;
    for (auto d : g.getDaughters()) {
      if (!d.isAvailable()) continue;
      if (isElectronPDG(d.getPDG())) ++nEle;
    }
    if (nEle < 2) return false;
    const V3 v = g.getVertex();
    const double r = std::hypot((double)v.x, (double)v.y);
    return (r >= rMinConv && r <= rMaxConv);
  }

  // Very generic 'reco object looks like a photon':
  // neutral, has cluster(s), no tracks.
  inline bool isPhotonLikeReco(const edm4hep::ReconstructedParticle& rp){
    if (std::abs(rp.getCharge()) > 1e-6) return false;
    if (rp.getClusters().empty())        return false;
    if (!rp.getTracks().empty())         return false;
    return true;
  }
}

struct PFOtoMCviaClusterLink final
  : k4FWCore::Consumer<void(const edm4hep::ReconstructedParticleCollection&,
                            const edm4hep::MCParticleCollection&,
                            const edm4hep::ClusterCollection&,
                            const edm4hep::RecoMCParticleLinkCollection&)> {
  PFOtoMCviaClusterLink(const std::string& name, ISvcLocator* svcLoc)
    : Consumer(name, svcLoc, {
        KeyValues("InputPFOs",        {"PandoraPFOs"}),
        KeyValues("InputMCParticles", {"MCParticles"}),
        KeyValues("InputClusters",    {"PandoraClusters"}),
        KeyValues("InputRecoMC",      {"MCTruthRecoLink"})
      }) {}

  // Histogram directory & binning
  Gaudi::Property<std::string> m_histPath{this,"HistPath","/MC/Cluster",
      "THistSvc directory (must start with the stream name, e.g. '/MC/...')."};
  Gaudi::Property<int>    m_binsE{this,"BinsE",64,"Bins in energy [GeV]."};
  Gaudi::Property<double> m_eMin{this,"EMin",0.,"Min E [GeV]."};
  Gaudi::Property<double> m_eMax{this,"EMax",128.,"Max E [GeV]."};
  Gaudi::Property<int>    m_binsCos{this,"BinsCosT",30,"Bins in cos(theta)."};

  // Converted-photon recognition knobs
  Gaudi::Property<float>  m_rMinConv{this,"ConvRmin",  5.0f,"Min R [mm] for conversion (photon vertex)."};
  Gaudi::Property<float>  m_rMaxConv{this,"ConvRmax",800.0f,"Max R [mm] for conversion (photon vertex)."};

  // Services
  mutable SmartIF<ITHistSvc> m_histSvc;
  mutable bool m_booked=false;

  // Generic spectra
  mutable TH1F *h_mcE=nullptr, *h_pfoE=nullptr, *h_mcCos=nullptr, *h_pfoCos=nullptr;
  // Response (Reco vs MC)
  mutable TH2F *h_respE=nullptr, *h_respCos=nullptr;
  // Efficiency as num/den (generic)
  mutable TH1F *h_effE_num=nullptr, *h_effE_den=nullptr;
  mutable TH1F *h_effCos_num=nullptr, *h_effCos_den=nullptr;
  // PID confusion (MC on X, Reco on Y)
  mutable TH2F *h_pidConf=nullptr;

  // --- NEW: photon split efficiencies (converted / unconverted)
  mutable TH1F *h_effE_den_g_conv=nullptr,  *h_effE_num_g_conv=nullptr;
  mutable TH1F *h_effE_den_g_unco=nullptr,  *h_effE_num_g_unco=nullptr;
  mutable TH1F *h_effC_den_g_conv=nullptr,  *h_effC_num_g_conv=nullptr;
  mutable TH1F *h_effC_den_g_unco=nullptr,  *h_effC_num_g_unco=nullptr;

  void bookIfNeeded() const {
    if (m_booked) return;
    m_histSvc = service("THistSvc");
    if (!m_histSvc) { error()<<"THistSvc not found"<<endmsg; return; }
    auto reg = [&](TH1* h){ return m_histSvc->regHist(m_histPath.value()+"/"+h->GetName(), h).isSuccess(); };

    // Generic
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

    h_pidConf    = new TH2F("pidConf","PID confusion;MC category;Reco category",
                             6, -0.5, 5.5, 6, -0.5, 5.5);

    // NEW: photon converted/unconverted splits (den/num vs E and vs cosθ)
    h_effE_den_g_conv  = new TH1F("effE_den_gamma_conv",  "Photon (converted) denominator vs E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_effE_num_g_conv  = new TH1F("effE_num_gamma_conv",  "Photon (converted) numerator vs E;E [GeV];Entries",   m_binsE, m_eMin, m_eMax);
    h_effE_den_g_unco  = new TH1F("effE_den_gamma_unconv","Photon (unconverted) denominator vs E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_effE_num_g_unco  = new TH1F("effE_num_gamma_unconv","Photon (unconverted) numerator vs E;E [GeV];Entries",   m_binsE, m_eMin, m_eMax);

    h_effC_den_g_conv  = new TH1F("effCos_den_gamma_conv",  "Photon (converted) denominator vs cos#theta;cos#theta;Entries", m_binsCos, -1, 1);
    h_effC_num_g_conv  = new TH1F("effCos_num_gamma_conv",  "Photon (converted) numerator vs cos#theta;cos#theta;Entries",   m_binsCos, -1, 1);
    h_effC_den_g_unco  = new TH1F("effCos_den_gamma_unconv","Photon (unconverted) denominator vs cos#theta;cos#theta;Entries", m_binsCos, -1, 1);
    h_effC_num_g_unco  = new TH1F("effCos_num_gamma_unconv","Photon (unconverted) numerator vs cos#theta;cos#theta;Entries",   m_binsCos, -1, 1);

    bool ok=true;
    ok&=reg(h_mcE); ok&=reg(h_pfoE); ok&=reg(h_mcCos); ok&=reg(h_pfoCos);
    ok&=reg(h_respE); ok&=reg(h_respCos);
    ok&=reg(h_effE_num); ok&=reg(h_effE_den); ok&=reg(h_effCos_num); ok&=reg(h_effCos_den);
    ok&=reg(h_pidConf);

    ok&=reg(h_effE_den_g_conv);  ok&=reg(h_effE_num_g_conv);
    ok&=reg(h_effE_den_g_unco);  ok&=reg(h_effE_num_g_unco);
    ok&=reg(h_effC_den_g_conv);  ok&=reg(h_effC_num_g_conv);
    ok&=reg(h_effC_den_g_unco);  ok&=reg(h_effC_num_g_unco);

    if (!ok) error()<<"Failed to register some histograms (check HistPath stream)."<<endmsg;

    m_booked=true;
  }

  void operator()(const edm4hep::ReconstructedParticleCollection& pfos,
                  const edm4hep::MCParticleCollection& mcps,
                  const edm4hep::ClusterCollection& /*clus*/,
                  const edm4hep::RecoMCParticleLinkCollection& links) const override {
    bookIfNeeded(); if (!m_booked) return;

    // MC index by ObjectID
    std::map<podio::ObjectID, size_t, ObjectIDLess> mcIndex;
    for (size_t i=0;i<mcps.size();++i) mcIndex[mcps[i].getObjectID()] = i;

    // Fill MC denominators/spectra once per event
    for (const auto& mc : mcps) {
      const float eMC  = mc.getEnergy();
      const float cMC  = safeCosTheta(mc.getMomentum());
      h_mcE->Fill(eMC);
      h_mcCos->Fill(cMC);
      h_effE_den->Fill(eMC);
      h_effCos_den->Fill(cMC);
    }

    // Group RecoMC links by reco id
    std::map<podio::ObjectID, std::vector<const edm4hep::RecoMCParticleLink*>, ObjectIDLess> byReco;
    byReco.clear();
    for (const auto& L : links) {
      if (!L.getFrom().isAvailable() || !L.getTo().isAvailable()) continue;  // from: Reco, to: MC
      byReco[L.getFrom().getObjectID()].push_back(&L);
    }

    // For filling photon split efficiencies we also want a "best MC per reco" and vice-versa decisions:
    // Here we select best linked MC per Reco by max weight, and best Reco per MC by max weight seen during the Reco loop.
    std::vector<int> bestRecoForMC(mcps.size(), -1);
    std::vector<float> bestWForMC(mcps.size(), -1.f);

    // Loop over PFOs: fill reco spectra, response, confusion & remember bestRecoForMC
    for (size_t iReco = 0; iReco < pfos.size(); ++iReco) {
      const auto& rp = pfos[iReco];

      const float eRec = rp.getEnergy();
      const float cRec = safeCosTheta(rp.getMomentum());
      h_pfoE->Fill(eRec);
      h_pfoCos->Fill(cRec);

      // pick best linked MC by weight
      int bestMCidx = -1; float bestW = -1.f;
      auto itL = byReco.find(rp.getObjectID());
      if (itL != byReco.end()) {
        for (auto Lptr : itL->second) {
          const auto mid = Lptr->getTo().getObjectID();
          auto itMC = mcIndex.find(mid);
          if (itMC == mcIndex.end()) continue;
          const float w = Lptr->getWeight();
          if (w > bestW) { bestW = w; bestMCidx = (int)itMC->second; }
        }
      }

      if (bestMCidx >= 0) {
        const auto& mc = mcps[bestMCidx];
        const float eMC = mc.getEnergy();
        const float cMC = safeCosTheta(mc.getMomentum());
        h_respE->Fill(eMC, eRec);
        h_respCos->Fill(cMC, cRec);

        // confusion matrix fill (categories by dictionary 0..5)
        const int mcCat = pidCategory(mc.getPDG());
        int recoCat = 5; // default other
        if (isPhotonLikeReco(rp)) recoCat = 1;
        else if (std::abs(rp.getCharge()) > 1e-6) recoCat = 4;
        h_pidConf->Fill(mcCat, recoCat);

        // remember best Reco index for this MC
        if (bestW > bestWForMC[bestMCidx]) {
          bestWForMC[bestMCidx] = bestW;
          bestRecoForMC[bestMCidx] = (int)iReco;
        }
      }
    }
    // --- Fill photon converted/unconverted efficiency splits (den/num) ---
    // Denominators: all MC photons (converted/unconverted) vs E and vs cosθ
    // Numerators: those whose best matching Reco looks like a photon (isPhotonLikeReco).
    for (size_t i=0;i<mcps.size();++i) {
      const auto& mc = mcps[i];
      if (mc.getPDG()!=22) continue; // photons only

      const float eMC = mc.getEnergy();
      const float cMC = safeCosTheta(mc.getMomentum());

      // vertex type is Vector3d in EDM4hep (use .x .y .z)
      const bool conv = isConvertedPhoton<edm4hep::Vector3d>(mc, m_rMinConv, m_rMaxConv);

      // denominators
      (conv ? h_effE_den_g_conv  : h_effE_den_g_unco)->Fill(eMC);
      (conv ? h_effC_den_g_conv  : h_effC_den_g_unco)->Fill(cMC);

      // numerators if best reco exists and is photon-like
      const int rIdx = bestRecoForMC[i];
      if (rIdx>=0 && rIdx<(int)pfos.size()) {
        const auto& rp = pfos[rIdx];
        if (isPhotonLikeReco(rp)) {
          (conv ? h_effE_num_g_conv  : h_effE_num_g_unco)->Fill(eMC);
          (conv ? h_effC_num_g_conv  : h_effC_num_g_unco)->Fill(cMC);
        }
      }
    }
  }
};

DECLARE_COMPONENT(PFOtoMCviaClusterLink)
