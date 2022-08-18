
cmsrel CMSSW_12_4_6

cd CMSSW_12_4_6/src/

cmsenv

git clone -b nanIssueInIso https://github.com/swagata87/EGMInfoDumper

scram b

cmsRun EGMInfoDumper/EGMGenericNtupler/python/ConfFile_cfg.py