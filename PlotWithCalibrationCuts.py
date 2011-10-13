import ROOT,os,sys
from array import array
ROOT.gSystem.Load("libEXOROOT")

def GetData(filename,hist):
  #input data
  t = ROOT.TChain("tree")
  t.Add(filename)
  nentries = t.GetEntries()

  print("nentries = " + str(nentries))
  for evtID in range(nentries):
    t.GetEntry(evtID)
    ED = t.EventBranch
    if (evtID%10000 == 0):
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
        hist.Fill(energy)

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
    hist = ROOT.TH1D("hist","hist",75,0,3500)
    GetData(sys.argv[1],hist)
    canvas = ROOT.TCanvas()
    canvas.cd()
    hist.Draw()
    canvas.Update()
    raw_input("enter to quit")
