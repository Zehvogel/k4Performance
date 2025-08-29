/*
 * PFOtoMCviaClusterLink: PFOs + Reco<->MC links, with 'selected' histograms
 * and runtime toggles for selection + debug counters.
 *
 * Selection (CLD-style, arXiv:1911.12230):
 *  - SelectByType:     reco/MC same coarse type (default true)
 *  - SelectByAngles:   |dphi| < MaxDeltaPhi && |dtheta| < MaxDeltaTheta  (defaults: 2 mrad / 1 mrad)
 *  - SelectPtForCharged: |pT(reco)/pT(true)-1| < MaxRelDpT (default 5%) for charged MC
 *
 * Debug:
 *  - DebugCounters: print rolling and final summaries
 *
 * Histograms (under HistPath, default "/PFO/PFOCluster"):
 *  Inclusive:  mcE, pfoE, mcCos, pfoCos, respE, effE_den/num, effCos_den/num, pidConf, respR
 *  Selected:   mcE_sel, pfoE_sel, mcCos_sel, pfoCos_sel, respE_sel,
 *              effE_sel_den/num, effCos_sel_den/num, respR_sel,
 *              dE_sel [GeV], dpT_sel [GeV] (charged MC only), dTheta_sel [rad], dPhi_sel [rad, signed]
 */

#include "k4FWCore/Consumer.h"

#include "Gaudi/Property.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/SmartIF.h"

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
#include <cstdint>

namespace {
  template <typename V3>
  inline float safeCosTheta(const V3& p) {
    const double px = (double)p.x, py = (double)p.y, pz = (double)p.z;
    const double mag = std::sqrt(px*px + py*py + pz*pz);
    return mag > 0. ? float(pz / mag) : 1.f;
  }
  template <typename V3>
  inline double thetaOf(const V3& p) {
    const double px = (double)p.x, py = (double)p.y, pz = (double)p.z;
    return std::atan2(std::sqrt(px*px + py*py), pz);
  }
  template <typename V3>
  inline double phiOf(const V3& p) { return std::atan2((double)p.y,(double)p.x); }
  template <typename V3>
  inline double pTOf(const V3& p) { return std::hypot((double)p.x,(double)p.y); }

  // absolute minimal deltaPhi in [0, pi]
  inline double deltaPhiAbs(double a, double b) {
    double d = std::fabs(a - b);
    if (d > M_PI) d = 2.0*M_PI - d;
    return d;
  }
  // signed minimal deltaPhi in (-pi, pi]
  inline double deltaPhiSigned(double a, double b) {
    double d = a - b;
    while (d <= -M_PI) d += 2.0*M_PI;
    while (d  >  M_PI) d -= 2.0*M_PI;
    return d;
  }

  struct ObjectIDLess {
    bool operator()(const podio::ObjectID& a, const podio::ObjectID& b) const noexcept {
      return (a.collectionID < b.collectionID) ||
             (a.collectionID == b.collectionID && a.index < b.index);
    }
  };
  // Coarse categories (align with Python): 0 mu, 1 gamma, 2 e, 3 nh, 4 ch, 5 other
  inline int pidBinMC(int pdg) {
    switch (pdg) {
      case 22: return 1;
      case 11: case -11: return 2;
      case 13: case -13: return 0;
      case 211: case -211: return 4;
      default: return 3;
    }
  }
  inline int pidBinReco(const edm4hep::ReconstructedParticle& rp) {
    if (std::fabs(rp.getCharge()) > 0.5) return 4; // charged hadron-like
    const double m = std::fabs((double)rp.getMass());
    if (m < 5e-3 && rp.getEnergy() > 0) return 1;  // gamma-like
    return 3; // other neutral
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

  // Output / binning
  Gaudi::Property<std::string> m_histPath{this,"HistPath","/PFO/PFOCluster",
      "THistSvc directory (must include stream, e.g. '/PFO/...')."};
  Gaudi::Property<int>    m_binsE {this,"BinsE",    64,   "Bins in energy [GeV]."};
  Gaudi::Property<double> m_eMin  {this,"EMin",     0.0,  "Min E [GeV]."};
  Gaudi::Property<double> m_eMax  {this,"EMax",     128., "Max E [GeV]."};
  Gaudi::Property<int>    m_binsC {this,"BinsCosT", 30,   "Bins in cos(theta)."};

  // Selection toggles & thresholds (these gate the 'selected' histograms)
  Gaudi::Property<bool>   m_DebugCounters     {this,"DebugCounters",      false, "Print rolling/final counters"};
  Gaudi::Property<bool>   m_SelectByType      {this,"SelectByType",       true,  "Require same coarse type"};
  Gaudi::Property<bool>   m_SelectByAngles    {this,"SelectByAngles",     true,  "Apply dphi/dtheta cuts"};
  Gaudi::Property<bool>   m_SelectPtForCharged{this,"SelectPtForCharged", true,  "Apply pT cut for charged"};
  Gaudi::Property<double> m_maxDphi           {this,"MaxDeltaPhi",        0.002, "Max |delta phi| in rad (2 mrad)"};
  Gaudi::Property<double> m_maxDtheta         {this,"MaxDeltaTheta",      0.001, "Max |delta theta| in rad (1 mrad)"};
  Gaudi::Property<double> m_maxRelDpT         {this,"MaxRelDpT",          0.05,  "Max |pT(reco)/pT(true)-1| for charged"};

  // Response ratio (Erec/Etrue - 1)
  Gaudi::Property<int>    m_binsRespR{this,"BinsRespR", 200, "Bins for response ratio"};
  Gaudi::Property<double> m_respRmin{this,"RespRMin", -0.5,  "Min response ratio"};
  Gaudi::Property<double> m_respRmax{this,"RespRMax", +0.5,  "Max response ratio"};

  // diff histogram ranges
  Gaudi::Property<int>    m_binsDiffE {this,"BinsDiffE",   200, "Bins for Delta E [GeV]"};
  Gaudi::Property<double> m_dEmin     {this,"DeltaEMin",   -20.,"Min Delta E [GeV]"};
  Gaudi::Property<double> m_dEmax     {this,"DeltaEMax",   +20.,"Max Delta E [GeV]"};

  Gaudi::Property<int>    m_binsDiffPt{this,"BinsDiffPt",  200, "Bins for Delta pT [GeV]"};
  Gaudi::Property<double> m_dPtMin    {this,"DeltaPtMin",  -20.,"Min Delta pT [GeV]"};
  Gaudi::Property<double> m_dPtMax    {this,"DeltaPtMax",  +20.,"Max Delta pT [GeV]"};

  Gaudi::Property<int>    m_binsDiffTh{this,"BinsDiffTh",  200, "Bins for Delta theta [rad]"};
  Gaudi::Property<double> m_dThMin    {this,"DeltaThMin",  -0.02,"Min Delta theta [rad]"};
  Gaudi::Property<double> m_dThMax    {this,"DeltaThMax",  +0.02,"Max Delta theta [rad]"};

  Gaudi::Property<int>    m_binsDiffPh{this,"BinsDiffPh",  200, "Bins for Delta phi [rad]"};
  Gaudi::Property<double> m_dPhMin    {this,"DeltaPhMin",  -0.02,"Min Delta phi [rad]"};
  Gaudi::Property<double> m_dPhMax    {this,"DeltaPhMax",  +0.02,"Max Delta phi [rad]"};

  // Services & booking flag
  mutable SmartIF<ITHistSvc> m_histSvc;
  mutable bool m_booked = false;

  // Inclusive
  mutable TH1F *h_mcE=nullptr, *h_pfoE=nullptr, *h_mcCos=nullptr, *h_pfoCos=nullptr;
  mutable TH2F *h_respE=nullptr;
  mutable TH1F *h_effE_den=nullptr, *h_effE_num=nullptr;
  mutable TH1F *h_effCos_den=nullptr, *h_effCos_num=nullptr;
  mutable TH2F *h_pidConf=nullptr;
  mutable TH1F *h_respR=nullptr;

  // Selected
  mutable TH1F *h_mcE_sel=nullptr, *h_pfoE_sel=nullptr;
  mutable TH1F *h_mcCos_sel=nullptr, *h_pfoCos_sel=nullptr;
  mutable TH2F *h_respE_sel=nullptr;
  mutable TH1F *h_effE_sel_den=nullptr, *h_effE_sel_num=nullptr;
  mutable TH1F *h_effCos_sel_den=nullptr, *h_effCos_sel_num=nullptr;
  mutable TH1F *h_respR_sel=nullptr;

  // Selected diffs
  mutable TH1F *h_dE_sel=nullptr, *h_dpT_sel=nullptr, *h_dTheta_sel=nullptr, *h_dPhi_sel=nullptr;

  // Debug counters
  mutable uint64_t m_evt{0}, m_nPrim{0}, m_nLinks{0}, m_nPrimWithReco{0}, m_nPrimSelected{0};

  void bookIfNeeded() const {
    if (m_booked) return;
    m_histSvc = service("THistSvc");
    if (!m_histSvc) { error() << "THistSvc not found" << endmsg; return; }

    auto reg = [&](TH1* h){
      const std::string full = m_histPath.value() + "/" + h->GetName();
      return m_histSvc->regHist(full, h).isSuccess();
    };

    // Inclusive
    h_mcE   = new TH1F("mcE","MC E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_pfoE  = new TH1F("pfoE","PFO E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_mcCos = new TH1F("mcCos","MC cos(theta);cos(theta);Entries", m_binsC, -1, 1);
    h_pfoCos= new TH1F("pfoCos","PFO cos(theta);cos(theta);Entries", m_binsC, -1, 1);

    h_respE = new TH2F("respE","Reco E vs MC E;MC E [GeV];Reco E [GeV]",
                       m_binsE, m_eMin, m_eMax, m_binsE, m_eMin, m_eMax);

    h_effE_den   = new TH1F("effE_den","Efficiency denominator vs E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_effE_num   = new TH1F("effE_num","Efficiency numerator vs E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_effCos_den = new TH1F("effCos_den","Efficiency denominator vs cos#theta;cos#theta;Entries", m_binsC, -1, 1);
    h_effCos_num = new TH1F("effCos_num","Efficiency numerator vs cos#theta;cos#theta;Entries", m_binsC, -1, 1);

    h_pidConf = new TH2F("pidConf","PID confusion;MC category;Reco category", 6, 0, 6, 6, 0, 6);
    // add labels to match the Python categories
    {
      static const char* labels[6] = {"mu","gamma","e","NH","CH","other"};
      for (int i = 1; i <= 6; ++i) {
        h_pidConf->GetXaxis()->SetBinLabel(i, labels[i-1]);
        h_pidConf->GetYaxis()->SetBinLabel(i, labels[i-1]);
      }
      h_pidConf->GetXaxis()->SetLabelSize(0.05);
      h_pidConf->GetYaxis()->SetLabelSize(0.05);
    }

    h_respR   = new TH1F("respR","Response (Erec/Etrue - 1);(Erec/Etrue - 1);Entries",
                         m_binsRespR, m_respRmin, m_respRmax);

    // Selected
    h_mcE_sel   = new TH1F("mcE_sel","MC E (selected);E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_pfoE_sel  = new TH1F("pfoE_sel","PFO E (selected);E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_mcCos_sel = new TH1F("mcCos_sel","MC cos(theta) (selected);cos(theta);Entries", m_binsC, -1, 1);
    h_pfoCos_sel= new TH1F("pfoCos_sel","PFO cos(theta) (selected);cos(theta);Entries", m_binsC, -1, 1);

    h_respE_sel = new TH2F("respE_sel","Reco E vs MC E (selected);MC E [GeV];Reco E [GeV]",
                           m_binsE, m_eMin, m_eMax, m_binsE, m_eMin, m_eMax);

    h_effE_sel_den   = new TH1F("effE_sel_den","Selected efficiency denominator vs E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_effE_sel_num   = new TH1F("effE_sel_num","Selected efficiency numerator vs E;E [GeV];Entries", m_binsE, m_eMin, m_eMax);
    h_effCos_sel_den = new TH1F("effCos_sel_den","Selected efficiency denominator vs cos#theta;cos#theta;Entries", m_binsC, -1, 1);
    h_effCos_sel_num = new TH1F("effCos_sel_num","Selected efficiency numerator vs cos#theta;cos#theta;Entries", m_binsC, -1, 1);

    h_respR_sel = new TH1F("respR_sel","Response (Erec/Etrue - 1) (selected);(Erec/Etrue - 1);Entries",
                           m_binsRespR, m_respRmin, m_respRmax);

    // Selected diffs
    h_dE_sel     = new TH1F("dE_sel",     "Delta E (reco - true) (selected);#Delta E [GeV];Entries",
                             m_binsDiffE, m_dEmin, m_dEmax);
    h_dpT_sel    = new TH1F("dpT_sel",    "Delta p_{T} (reco - true) (selected, charged MC only);#Delta p_{T} [GeV];Entries",
                             m_binsDiffPt, m_dPtMin, m_dPtMax);
    h_dTheta_sel = new TH1F("dTheta_sel", "Delta #theta (reco - true) (selected);#Delta #theta [rad];Entries",
                             m_binsDiffTh, m_dThMin, m_dThMax);
    h_dPhi_sel   = new TH1F("dPhi_sel",   "Delta #phi (signed) (selected);#Delta #phi [rad];Entries",
                             m_binsDiffPh, m_dPhMin, m_dPhMax);

    bool ok = true;
    ok&=reg(h_mcE); ok&=reg(h_pfoE); ok&=reg(h_mcCos); ok&=reg(h_pfoCos);
    ok&=reg(h_respE);
    ok&=reg(h_effE_den); ok&=reg(h_effE_num); ok&=reg(h_effCos_den); ok&=reg(h_effCos_num);
    ok&=reg(h_pidConf); ok&=reg(h_respR);
    ok&=reg(h_mcE_sel); ok&=reg(h_pfoE_sel); ok&=reg(h_mcCos_sel); ok&=reg(h_pfoCos_sel);
    ok&=reg(h_respE_sel);
    ok&=reg(h_effE_sel_den); ok&=reg(h_effE_sel_num); ok&=reg(h_effCos_sel_den); ok&=reg(h_effCos_sel_num);
    ok&=reg(h_respR_sel);

    // register diffs
    ok&=reg(h_dE_sel); ok&=reg(h_dpT_sel); ok&=reg(h_dTheta_sel); ok&=reg(h_dPhi_sel);

    if (!ok) error() << "Failed to register some histograms (check HistPath)" << endmsg;
    m_booked = true;
  }

  // Predicates with toggles
  bool sameType(const edm4hep::ReconstructedParticle& rp, const edm4hep::MCParticle& mc) const {
    if (!m_SelectByType) return true;
    return pidBinReco(rp) == pidBinMC(mc.getPDG());
  }
  bool passAngles(const edm4hep::ReconstructedParticle& rp, const edm4hep::MCParticle& mc) const {
    if (!m_SelectByAngles) return true;
    const double dphi = deltaPhiAbs(phiOf(rp.getMomentum()), phiOf(mc.getMomentum())); // radians
    const double dth  = std::fabs(thetaOf(rp.getMomentum()) - thetaOf(mc.getMomentum())); // radians
    return (dphi < m_maxDphi) && (dth < m_maxDtheta);
  }
  bool passPtIfCharged(const edm4hep::ReconstructedParticle& rp, const edm4hep::MCParticle& mc) const {
    if (!m_SelectPtForCharged) return true;
    if (std::fabs(mc.getCharge()) < 0.5) return true; // only applies for charged MC
    const double prt = pTOf(rp.getMomentum());
    const double pmt = pTOf(mc.getMomentum());
    if (pmt <= 0) return false;
    return std::fabs(prt / pmt - 1.0) < m_maxRelDpT;
  }

  void operator()(const edm4hep::ReconstructedParticleCollection& pfos,
                  const edm4hep::MCParticleCollection& mcps,
                  const edm4hep::ClusterCollection& /*clus*/,
                  const edm4hep::RecoMCParticleLinkCollection& links) const override {

    bookIfNeeded();
    if (!m_booked) return;

    // Index MC by ObjectID
    std::map<podio::ObjectID, size_t, ObjectIDLess> mcIndex;
    for (size_t i=0; i<mcps.size(); ++i) mcIndex[ mcps[i].getObjectID() ] = i;

    // Gun primaries (fallback to all)
    std::vector<size_t> primaries;
    primaries.reserve(mcps.size());
    for (size_t i=0; i<mcps.size(); ++i) if (mcps[i].parents_size()==0) primaries.push_back(i);
    if (primaries.empty()) for (size_t i=0;i<mcps.size();++i) primaries.push_back(i);

    // Inclusive truth spectra
    for (const auto& mc : mcps) {
      const float eMC = mc.getEnergy();
      const float cMC = safeCosTheta(mc.getMomentum());
      h_mcE->Fill(eMC); h_mcCos->Fill(cMC);
    }

    // Selected denominators once per primary (efficiency baseline)
    for (auto idx : primaries) {
      const auto& mc = mcps[idx];
      if (h_effE_sel_den)   h_effE_sel_den->Fill(mc.getEnergy());
      if (h_effCos_sel_den) h_effCos_sel_den->Fill(safeCosTheta(mc.getMomentum()));
      m_nPrim++;
    }

    // Group links by Reco (from=Reco, to=MC)
    std::map<podio::ObjectID, std::vector<const edm4hep::RecoMCParticleLink*>, ObjectIDLess> linksByReco;
    for (const auto& L : links) {
      if (!L.getFrom().isAvailable() || !L.getTo().isAvailable()) continue;
      linksByReco[L.getFrom().getObjectID()].push_back(&L);
    }

    // Inclusive reco spectra and inclusive response + inclusive efficiency (numerator = link exists)
    for (size_t ir=0; ir<pfos.size(); ++ir) {
      const auto& rp = pfos[ir];
      const float eRec = rp.getEnergy();
      const float cRec = safeCosTheta(rp.getMomentum());
      h_pfoE->Fill(eRec); h_pfoCos->Fill(cRec);

      auto itL = linksByReco.find(rp.getObjectID());
      if (itL != linksByReco.end()) {
        int bestIdx = -1; float bestW = -1.f;
        for (auto Lptr : itL->second) {
          auto itMC = mcIndex.find(Lptr->getTo().getObjectID());
          if (itMC == mcIndex.end()) continue;
          const float w = Lptr->getWeight();
          if (w > bestW) { bestW = w; bestIdx = (int)itMC->second; }
          m_nLinks++;
        }
        if (bestIdx >= 0) {
          const auto& mc = mcps[bestIdx];
          const float eMC = mc.getEnergy();
          const float cMC = safeCosTheta(mc.getMomentum());
          if (h_respE)   h_respE->Fill(eMC, eRec);
          if (h_effE_num)   h_effE_num->Fill(eMC);   // numerator counts existence of link
          if (h_effCos_num) h_effCos_num->Fill(cMC); // numerator counts existence of link
          if (h_respR && eMC > 0.f) h_respR->Fill((eRec / eMC) - 1.f);
          if (h_pidConf) h_pidConf->Fill(pidBinMC(mc.getPDG()), pidBinReco(rp));
        }
      }
    }

    // Best reco per MC (by weight)
    std::vector<int>   bestRecoForMC(mcps.size(), -1);
    std::vector<float> bestWforMC(mcps.size(), -1.f);
    for (size_t ir = 0; ir < pfos.size(); ++ir) {
      auto itL = linksByReco.find(pfos[ir].getObjectID());
      if (itL == linksByReco.end()) continue;
      for (auto Lptr : itL->second) {
        auto itMC = mcIndex.find(Lptr->getTo().getObjectID());
        if (itMC == mcIndex.end()) continue;
        const size_t imc = itMC->second;
        const float w = Lptr->getWeight();
        if (w > bestWforMC[imc]) { bestWforMC[imc] = w; bestRecoForMC[imc] = (int)ir; }
      }
    }

    // Selected fills once per primary (apply toggles here)
    m_evt++;
    for (auto imc : primaries) {
      int ir = bestRecoForMC[imc];
      if (ir < 0) continue; // no reco for this primary
      m_nPrimWithReco++;

      const auto& mc = mcps[imc];
      const auto& rp = pfos[(size_t)ir];

      const float eMC  = mc.getEnergy();
      const float cMC  = safeCosTheta(mc.getMomentum());
      const float eRec = rp.getEnergy();
      const float cRec = safeCosTheta(rp.getMomentum());

      const bool okType = sameType(rp, mc);
      const bool okAng  = passAngles(rp, mc);         // uses radians
      const bool okPt   = passPtIfCharged(rp, mc);    // DeltapT relative for charged MC

      if (okType && okAng && okPt) {
        m_nPrimSelected++;

        if (h_effE_sel_num)   h_effE_sel_num->Fill(eMC);
        if (h_effCos_sel_num) h_effCos_sel_num->Fill(cMC);

        if (h_respE_sel)   h_respE_sel->Fill(eMC, eRec);
        if (h_respR_sel && eMC > 0.f) h_respR_sel->Fill((eRec / eMC) - 1.f);

        if (h_mcE_sel)   h_mcE_sel->Fill(eMC);
        if (h_pfoE_sel)  h_pfoE_sel->Fill(eRec);
        if (h_mcCos_sel) h_mcCos_sel->Fill(cMC);
        if (h_pfoCos_sel)h_pfoCos_sel->Fill(cRec);

        // Diffs for selected pairs
        const double thReco = thetaOf(rp.getMomentum());
        const double thMC   = thetaOf(mc.getMomentum());
        const double phReco = phiOf(rp.getMomentum());
        const double phMC   = phiOf(mc.getMomentum());
        const double dE     = (double)eRec - (double)eMC;
        const double dTh    = thReco - thMC;                    // radians
        const double dPh    = deltaPhiSigned(phReco, phMC);     // radians, signed, circular
        if (h_dE_sel)     h_dE_sel->Fill(dE);
        if (h_dTheta_sel) h_dTheta_sel->Fill(dTh);
        if (h_dPhi_sel)   h_dPhi_sel->Fill(dPh);

        // Only for charged MC: fill ΔpT
        if (std::fabs(mc.getCharge()) > 0.5) {
          const double dpt = pTOf(rp.getMomentum()) - pTOf(mc.getMomentum());
          if (h_dpT_sel) h_dpT_sel->Fill(dpt);
        }
      }
    }

    if (m_DebugCounters && (m_evt % 100 == 0)) {
      info() << "evt=" << m_evt
             << " primaries_so_far=" << m_nPrim
             << " primWithReco_so_far=" << m_nPrimWithReco
             << " primSelected_so_far=" << m_nPrimSelected
             << " links_so_far=" << m_nLinks
             << endmsg;
    }
  }

  StatusCode finalize() override {
    if (m_DebugCounters) {
      info() << "[Summary] events=" << m_evt
             << " primaries=" << m_nPrim
             << " primWithReco=" << m_nPrimWithReco
             << " primSelected=" << m_nPrimSelected
             << " links=" << m_nLinks
             << endmsg;
    }
    return Consumer::finalize();
  }
};

DECLARE_COMPONENT(PFOtoMCviaClusterLink)
