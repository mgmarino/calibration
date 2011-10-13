import ROOT,os,sys
from array import array
ROOT.gSystem.Load("libEXOROOT")
ROOT.gSystem.Load("fit2.so")

def GetData(filename):
  #input data
  t = ROOT.TChain("tree")
  t.Add(filename)
  nentries = t.GetEntries()

  print("nentries = " + str(nentries))
  for evtID in range(nentries):
    t.GetEntry(evtID)
    ED = t.EventBranch
    if (evtID%1000 == 0):
      print(str(evtID) +  " events processed ")

    nsc = ED.GetNumScintillationClusters()
    for scID in range(nsc):
      scint_cluster = ED.GetScintillationCluster(scID)
      if (scint_cluster.fTime > 1928000):
        continue 
      ncl = scint_cluster.GetNumChargeClusters()
      if(ncl != 1): # only single site events
        continue

      energy = 0.0
      good = True
      for clID in range(ncl):
        cc = scint_cluster.GetChargeClusterAt(clID);
        energy += cc.fPurityCorrectedEnergy
        good = IsFiducial(cc.fX,cc.fY,cc.fZ)
        if(not good):
          break
      if(good and energy > 2000. and energy < 3500.):
        #print("adding energy " + str(energy))
        ROOT.FIT.AddDataPoint(energy)

def IsFiducial(x,y,z):
  fiducial = True
  if (ROOT.TMath.Sqrt(x**2 + y**2) > 163):
    fiducial = False
  if (z > 172 or z < -172):
    fiducial = False
  if (z > -20 and z < 20):
    fiducial = False
  return fiducial

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("wrong number of arguments")
  else:
    ROOT.FIT.SetEmin(2000.0)
    ROOT.FIT.SetEmax(3500.0)
    GetData(sys.argv[1])
    params = array("d",[3.0,2700.0,100.0])
    ROOT.FIT.SetParams(params)
    canvas1 = ROOT.TCanvas()
    canvas2 = ROOT.TCanvas()
    canvas1.cd()
    hist = ROOT.FIT.GetHist("InitialGuess")
    hist.Draw()
    canvas1.Update()
    raw_input("enter to continue")
    del hist
    ROOT.FIT.fit()
    canvas2.cd()
    hist = ROOT.FIT.GetHist("Fit")
    hist.Draw()
    canvas2.Update()
    raw_input("enter to quit")
