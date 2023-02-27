# EGMInfoDumper

Generic EGM ntuplizer for RunIII, for multi purpose use 

```
cmsrel CMSSW_12_6_4
cd CMSSW_12_6_4/src/
cmsenv
git cms-init  
git clone -b run3ID git@github.com:swagata87/EgammaPostRecoTools.git  EgammaUser/EgammaPostRecoTools
git clone -b example_run3 https://github.com/swagata87/EGMInfoDumper
scram b
cmsRun EGMInfoDumper/EGMGenericNtupler/python/ConfFile_cfg.py
```

From your analysis setup you need to run Egamma post-reco tool like this:
```
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       era='run3',
                       runEnergyCorrections=False,
                       runVID=True,
)  

#Run Egamma post reco tool
process.p = cms.Path(process.egammaPostRecoSeq * process.demo)
```