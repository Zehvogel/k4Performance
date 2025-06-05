from Gaudi.Configuration import INFO, DEBUG

# Loading some files with MCParticles and clusters for testing purposes
from k4FWCore import IOSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = "/afs/cern.ch/user/b/brfranco/work/public/FCC-config/FCCee/FullSim/ALLEGRO/ALLEGRO_o1_v03/ALLEGRO_sim_digi_reco.root"
io_svc.Output = "ALLEGRO_sim_digi_reco_perf.root"

from Configurables import CaloClusterMCParticleLinker
CaloClusterMCParticleLinker = CaloClusterMCParticleLinker("CaloClusterMCParticleLinker",
        InputMCParticles = ["MCParticles"],
        InputCaloClusters = ["CalibratedEMBCaloClusters"],
        OutputCaloClusterMCParticleLink = ["CalibratedEMBCaloClusters_MCParticles_Link"],
        ProduceValidationPlots = True,
        cutoff_dR = 0.3,
        )


# Set auditor service
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
CaloClusterMCParticleLinker.AuditExecute = True

from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr
ApplicationMgr(
    TopAlg=[CaloClusterMCParticleLinker],
    EvtSel='NONE',
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
)
