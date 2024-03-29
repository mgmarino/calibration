import ROOT,os,sys
from array import array
from Utilities import FiducialCut

ROOT.gSystem.Load("libEXOROOT")
ROOT.gSystem.Load("fit2.so")

class DataProcessor:
  def __init__(self,min,max):
    self.calib = ROOT.EXOCalibManager.GetCalibManager()
    self.calib.SetMetadataAccessType("mysql")
    self.calib.SetPort(3306)
    self.calib.SetUser("rd_exo_cond_ro")
    self.calib.SetHost("mysql-node03.slac.stanford.edu")
    self.rec = ROOT.EXOReconstructionModule()
    self.wiregain = ROOT.EXOWireGainModule()
    self.pur = ROOT.EXOLifetimeCalibModule()
    self.rec.set_use_true_t0_flag(False)
    self.rec.set_ATeam_YMatch_flag(1)
    self.rec.set_pattern_recognition(4)
    self.rec.set_SumBothAPDPlanes(False)
    #self.rec.set_ChannelFitChiSquareCut(-999.9)
    self.pur.SetFlavor("vanilla")
    self.Beginofrun = True
    self.Emin = min
    self.Emax = max
    self.rec.Initialize()
    self.wiregain.Initialize()
    self.pur.Initialize()
    self.Cut = FiducialCut()

  def ProcessEvent(self,ED):
    ED.GetWaveformData().Decompress()
    if(self.Beginofrun):
      self.Beginofrun = False
      self.rec.BeginOfRun(ED)
      self.wiregain.BeginOfRun(ED)
      self.pur.BeginOfRun(ED)
    self.rec.ProcessEvent(ED)
    self.wiregain.ProcessEvent(ED)
    self.pur.ProcessEvent(ED)

  def SetLowerWindow(self,win):
    self.rec.rec.set_lower_fit_boundary_wire_microseconds(win)

  def SetUpperWindow(self,win):
    self.rec.rec.set_upper_fit_boundary_wire_microseconds(win)

  def AddEnergiesAfterCut(self,ED,multisite = False):
    nscl = ED.GetNumScintillationClusters()
    for i in range(nscl):
      sc = ED.GetScintillationCluster(i)
      nccl = sc.GetNumChargeClusters()
      if (nccl != 1) and not multisite:
        continue
      if (nccl <= 1) and multisite:
        continue
      energy = 0.0
      good = True
      for j in range(nccl):
        cc = sc.GetChargeClusterAt(j)
        good = self.Cut.IsFiducial(cc.fX,cc.fY,cc.fZ)
        if not good:
          break
        energy += cc.fPurityCorrectedEnergy
        #energy += cc.fCorrectedEnergy
      if good and (energy > self.Emin) and (energy < self.Emax):
        ROOT.FIT.AddDataPoint(energy)

def ProcessRun(runnumber,lowerWin,upperWin,Emin,Emax,initialParams,multisite):
  ROOT.FIT.ClearData()
  params = array("d",initialParams)
  ROOT.FIT.SetParams(params)
  processor = DataProcessor(Emin,Emax)
  processor.SetLowerWindow(lowerWin)
  processor.SetUpperWindow(upperWin)
  basedir = "$EXOROOTDATA/"
  path = os.path.expandvars(basedir + str(runnumber))
  t = ROOT.TChain("tree")
  #path += "/run*.root"
  path += "/run00001926-000.root"
  #print(path)
  t.Add(path)
  nentries = t.GetEntries()
  print("Processing " + str(nentries) + " events")
  for i in range(nentries):
    if(not i%1000):
      print("event " + str(i))
    t.GetEntry(i)
    ED = t.EventBranch
    processor.ProcessEvent(ED)
    processor.AddEnergiesAfterCut(ED,multisite)
  ROOT.FIT.fit()
  params = ROOT.FIT.GetParams()
  errors = ROOT.FIT.GetParamErrors()
  return params[0],params[1]

def main(runs):
  print runs
  for run in runs:
    windows = [-20,-10,0,10,20,30]
    energies = []
    sigmas = []
    for window in windows:
      InitialParams = [3,2700,100]
      e, s =  ProcessRun(run,window,90,2000,3500,InitialParams,False)
      energies.append(e)
      sigmas.append(s)
    graph = ROOT.TGraphErrors(len(windows),array("d",windows),array("d",energies),array("d",len(windows)*[0]),array("d",sigmas))
    graph.Draw("AP")
    raw_input("hit enter to continue")


if __name__ == "__main__":
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    print("please specify run number(s)")
