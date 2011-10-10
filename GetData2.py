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

      # safe all charge clusters in array
      xcl = [] 
      ycl = [] 
      zcl = [] 
      ucl = [] 
      vcl = [] 
      dtcl = [] 
      epcl = [] 
      tcl = []
      dhalfcl = []
      nUWires = []

      for clID in range(ncl):
        charge_cluster = scint_cluster.GetChargeClusterAt(clID);

        xcl.append(charge_cluster.fX)
        ycl.append(charge_cluster.fY)
        zcl.append(charge_cluster.fZ)
        ucl.append(charge_cluster.fU)
        vcl.append(charge_cluster.fV)
        dtcl.append(charge_cluster.fDriftTime / 1000.0)
        epcl.append(charge_cluster.fPurityCorrectedEnergy)
        tcl.append(charge_cluster.fCollectionTime)
        dhalfcl.append(charge_cluster.fDetectorHalf)

        nU = charge_cluster.GetNumUWireSignals()
        nUWires.append(nU)

      #apply fiducial cut
      fiducial = True
      energy = 0.0
      for i in range(ncl):
        energy += epcl[i]
        if (ROOT.TMath.Sqrt(xcl[i]*xcl[i] + ycl[i]*ycl[i]) > 163):
          fiducial = False
          break
        if (zcl[i] > 172 or zcl[i] < -172):
          fiducial = False
          break
        if (zcl[i] > -20 and zcl[i] < 20):
          fiducial = False
          break
        #if (epcl[i] < 100.0):    # for Cs137 every charge cluster must at least have 100 keV
        #  fiducial = False
        #  continue 
      if (ncl == 0):
        fiducial = False
        continue
      if (ncl != 1):
        continue
      if (not fiducial):
        continue
      if (ncl == 0):
        continue
      if (energy < 2000.0 or energy > 3500.0):
        continue
      ROOT.AddDataPoint(energy)

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("wrong number of arguments")
  else:
    GetData(sys.argv[1])
    ROOT.SetEmin(2000.0)
    ROOT.SetEmax(3500.0)
    params = array("d",[3.8,2750.0,100.0])
    ROOT.SetParams(params)
    canvas = ROOT.TCanvas()
    canvas.cd()
    hist = ROOT.GetHist()
    hist.Draw()
    canvas.Update()
    raw_input("enter to continue")
    ROOT.fit()
    ROOT.GetHist()
    hist = ROOT.GetHist()
    hist.Draw()
    canvas.Update()
    raw_input("enter to quit")
