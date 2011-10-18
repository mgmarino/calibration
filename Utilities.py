import math,ROOT

class FiducialCut:
  def __init__(self,r=163.,lz=20.,uz=172.):
    self.rcut = r
    self.lowerzcut = lz
    self.upperzcut = uz

  def IsFiducial(x,y,z)
    if math.hypot(x,y) > self.rcut:
      return False
    if abs(z) > self.upperzcut:
      return False
    if abs(z) < self.lowerzcut:
      return False
    return True

class Reconstruction:
  def __init__(self):
    self.rec = ROOT.EXOReconstructionModule()
    self.rec.set_use_true_t0_flag(False)
    self.rec.set_ATeam_YMatch_flag(1)
    self.rec.set_pattern_recognition(4)
    self.rec.set_SumBothAPDPlanes(False)
    self.rec.Initialize()
    self.Beginofrun = True

  def ProcessEvent(self,ED):
    ED.GetWaveformData().Decompress()
    if(self.Beginofrun):
      self.Beginofrun = False
      self.rec.BeginOfRun(ED)
    self.rec.ProcessEvent(ED)

