# EGMInfoDumper

Generic EGM ntuplizer for RunII or RunIII, for multi purpose use 

```
cmsrel CMSSW_12_6_4
cd CMSSW_12_6_4/src/
cmsenv
git cms-init  
#Get the post reco tool, 
#if you need Run3 ID then get it the following way
git clone -b run3ID git@github.com:swagata87/EgammaPostRecoTools.git  EgammaUser/EgammaPostRecoTools
#Or, for Run2 UL, instead follow this https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
#Get example ntuplizer
git clone -b aod https://github.com/swagata87/EGMInfoDumper
scram b
cmsRun EGMInfoDumper/EGMGenericNtupler/python/ConfFile_cfg.py
```

From your analysis setup you need to run Egamma post-reco tool like this:
```
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       era='2018-UL',
                       runVID=True,
		       isMiniAOD=False
)  

#Run Egamma post reco tool
process.p = cms.Path(process.egammaPostRecoSeq * process.demo)
```