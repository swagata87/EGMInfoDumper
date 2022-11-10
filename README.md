
cmsrel CMSSW_12_4_6

cd CMSSW_12_4_6/src/

cmsenv

git-cms-init

git-cms-addpkg RecoEgamma/PhotonIdentification

git clone -b ExampleOfRhoCorrection https://github.com/swagata87/EGMInfoDumper

scram b

cmsRun EGMInfoDumper/EGMGenericNtupler/python/ConfFile_cfg.py