import ROOT,os,sys
from array import array

ROOT.gSystem.Load("libEXOROOT")
ROOT.gSystem.Load("fit2.so")

class DataProcessor:
  def __init__(self,min,max):
    self.rec = ROOT.EXOReconstructionModule()
    self.wiregain = ROOT.EXOWireGainModule()
    self.pur = ROOT.EXOLifetimeCalibModule()
    self.Beginofrun = True
    self.Emin = min
    self.Emax = max
    self.rec.Initialize()
    self.rec.rec.set_use_trigger_t0(False)
    self.rec.rec.set_ATeam_YMatch(1)
    self.rec.rec.set_pattern_recognition_apd(4)
    self.rec.rec.set_pattern_recognition_uwire(4)
    self.rec.rec.SumBothAPDPlanes = False
    self.rec.rec.ChannelFitChiSquareCut = -999.9
    self.wiregain.Initialize()
    self.pur.Initialize()
    self.pur.SetFlavor("vanilla")

  def ProcessEvent(self,ED):
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
        good = IsFiducial(cc.fX,cc.fY,cc.fZ)
        if not good:
          break
        energy += cc.fPurityCorrectedEnergy
      if good and (energy > self.Emin) and (energy < self.Emax):
        ROOT.FIT.AddDataPoint(energy)

def ProcessRun(runnumber,lowerWin,upperWin,Emin,Emax,initialParams,multisite):
  ROOT.FIT.ClearData()
  params = array("d",initialParams)
  ROOT.FIT.SetParams(params)
  processor = DataProcessor(Emin,Emax)
  processor.SetLowerWindow(lowerWin)
  processor.SetUpperWindow(upperWin)
  path = os.path.expandvars("$EXODATA/WIPP/root/"+str(runnumber))
  t = ROOT.TChain("tree")
  t.Add(path + "/run*.root")
  nentries = t.GetEntries()
  print("Processing " + str(nentries) + " events")
  for i in range(nentries):
    print("event " + str(i))
    t.GetEntry(i)
    ED = t.EventBranch
    processor.ProcessEvent(ED)
    processor.AddEnergiesAfterCut(ED,multisite)
  ROOT.FIT.fit()
  params = ROOT.FIT.GetParams()
  errors = ROOT.FIT.GetParamErrors()
  return params[0],params[1]

def IsFiducial(x,y,z):
  fiducial = True
  if (ROOT.TMath.Sqrt(x**2 + y**2) > 163):
    fiducial = False
  if (z > 172 or z < -172):
    fiducial = False
  if (z > -20 and z < 20):
    fiducial = False
  return fiducial

def main():
  runs = [1926]
  #runs = [2448]
  for run in runs:
    windows = [-20,-10,0,10,20,30]
    energies = []
    sigmas = []
    for window in windows:
      InitialParams = [3,2700,100]
      e, s =  ProcessRun(run,window,90,2000,3500,InitialParams,False)
      energies.append(e)
      sigmas.append(s)
    graph = ROOT.TGraphErrors(windows,energies,len(windows)*[0],sigmas)
    graph.Draw("AP")
    raw_input("hit enter to continue")


if __name__ == "__main__":
  main()
