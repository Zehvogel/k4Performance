#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from Gaudi.Configuration import DEBUG
from k4FWCore import ApplicationMgr
from k4FWCore import IOSvc
from Configurables import MCvalidation, THistSvc, EventDataSvc
from Configurables import PFOtoMCviaClusterLink, PFOtoMCviaCellLink
from k4FWCore.parseArgs import parser
parser_group = parser.add_argument_group("runPFOanalysis.py custom options")
parser_group.add_argument("-i","--inputFiles", action="extend", nargs="+", metavar=("file1", "file2"), help="One or multiple input files")
parser_group.add_argument("-o","--outputBasename", help="Basename of the output file(s)", default="PFO_analysis_output")
parsed_args = parser.parse_known_args()[0]

iosvc = IOSvc("IOSvc")
iosvc.Input = parsed_args.inputFiles

iosvc.CollectionNames = ["MCParticles", "PandoraPFOs", "PandoraClusters", "MCTruthRecoLink"]

clu = PFOtoMCviaClusterLink("ViaCluster")
clu.OutputLevel = DEBUG
clu.HistPath = "/PFO/PFOCluster"
clu.InputPFOs     = ["PandoraPFOs"]
clu.InputMCParticles = ["MCParticles"]
clu.InputClusters = ["PandoraClusters"]
clu.InputRecoMC   = ["MCTruthRecoLink"]
clu.DebugCounters = True
clu.SelectByType = False
clu.SelectByAngles = False
clu.SelectPtForCharged = False
# thresholds (radians/fraction)
#clu.MaxDeltaPhi = 0.002
#clu.MaxDeltaTheta = 0.001
#clu.MaxRelDpT = 0.05


# Histogram output (THistSvc)
THistSvc().Output = [f"PFO DATAFILE='{parsed_args.outputBasename}.root' OPT='RECREATE' TYP='ROOT'"]
ApplicationMgr(TopAlg=[clu],
               EvtSel="NONE",
               EvtMax=-1,
               ExtSvc=[EventDataSvc("EventDataSvc"), THistSvc("THistSvc"),iosvc],
               OutputLevel=DEBUG,
               )
