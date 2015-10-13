import datetime

beginYear = 1981
endYear = 1981
forcingInput = "/tigress/nwanders/Scripts/Seasonal/"
modelS = ["CanCM3","CanCM4", "FLOR"]
precVarNameS = ["prlr","prlr","pr"]
tempVarName = "tas"
refInput = "/tigress/nwanders/Scripts/Seasonal/refData/"
refModelS = ["PGF"] #["PGF", "CFS"]
precRefVarName = "prec"
tempRefVarName = "tas"
precCorFactor = 1.0
pctlInput = "/tigress/nwanders/Scripts/Seasonal/resultsNetCDF/"

master = open("/tigress/nwanders/Scripts/hydroSeasonal/jobs/VICmaster.sh", "w")

for m in range(len(modelS)):
  model = modelS[m]
  precVarName = precVarNameS[m]
  for r in range(len(refModelS)):
    print r
    refModel = refModelS[r]
    for year in range(beginYear,endYear+1):
      for month in range(1,13):
        startTime = '%d%.2d01' %(year, month)
        endTimeTemp = datetime.datetime(year+1, month, 1) - datetime.timedelta(days=1)
        endTime = '%d%.2d%.2d' %(endTimeTemp.year, endTimeTemp.month, endTimeTemp.day)
 
        job = open("/tigress/nwanders/Scripts/hydroSeasonal/jobs/VIC_"+model+"_"+refModel+"_"+startTime+".sh", "w")
        job.writelines("#!/bin/bash\n")
        job.writelines("#SBATCH -n 1   # node count\n")
        job.writelines("#SBATCH -t 23:59:59\n")
        job.writelines("#SBATCH --mail-type=fail\n")
        job.writelines("#SBATCH --mail-user=nwanders@princeton.edu\n")
        job.writelines("cd /tigress/nwanders/Scripts/hydroSeasonal/\n")
        job.writelines("python /tigress/nwanders/Scripts/hydroSeasonal/seasonScripts/runSeasonalVIC.py %s %s %s %s" %(startTime, endTime, model, refModel))
        job.close()
    
        master.writelines("sbatch /tigress/nwanders/Scripts/hydroSeasonal/jobs/VIC_"+model+"_"+refModel+"_"+startTime+".sh\n")

master.close()
