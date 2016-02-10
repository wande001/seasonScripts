start = "#!/bin/bash
# parallel job using 48 cores. and runs for 4 hours (max)
#SBATCH -n 1   # node count
#SBATCH -t 23:59:59
# sends mail when process begins, and
# when it ends. Make sure you define your email
# address.
#SBATCH --mail-type=fail
#SBATCH --mail-user=nwanders@princeton.edu

cd /tigress/nwanders/Scripts/hydroSeasonal/seasonScripts\n"

runs = list()

#runs[[1]] ="python calcMetricDischargeVIC.py 0 1 24 1 $lag $model $ref $varName"
#runs[[2]] = "python calcMetricDischargeVIC.py 1 2 24 2 $lag $model $ref $varName"
#runs[[3]] = "python calcMetricDischargeVIC.py 2 2 24 4 $lag $model $ref $varName"
#runs[[4]] = "python calcMetricDischargeVIC.py 3 2 24 6 $lag $model $ref $varName"
#runs[[5]] = "python calcMetricDischargeVIC.py 4 2 24 8 $lag $model $ref $varName"
#runs[[6]] = "python calcMetricDischargeVIC.py 5 2 24 10 $lag $model $ref $varName"
#runs[[7]] = "python calcMetricDischargeVIC.py 6 2 24 12 $lag $model $ref $varName"
#runs[[8]] = "python calcMetricDischargeVIC.py 7 2 24 14 $lag $model $ref $varName"
#runs[[9]] = "python calcMetricDischargeVIC.py 8 2 24 16 $lag $model $ref $varName"
#runs[[10]] = "python calcMetricDischargeVIC.py 9 2 24 18 $lag $model $ref $varName"
#runs[[11]] = "python calcMetricDischargeVIC.py 10 2 24 20 $lag $model $ref $varName"
#runs[[12]] = "python calcMetricDischargeVIC.py 11 2 24 22 $lag $model $ref $varName"
#runs[[13]] = "python calcMetricDischargeVIC.py 12 2 24 24 $lag $model $ref $varName"

runs[[1]] ="python calcMetricCatchments.py 0 1 24 1 $lag $model $ref $varName"
runs[[2]] = "python calcMetricCatchments.py 1 2 24 2 $lag $model $ref $varName"
runs[[3]] = "python calcMetricCatchments.py 2 2 24 4 $lag $model $ref $varName"
runs[[4]] = "python calcMetricCatchments.py 3 2 24 6 $lag $model $ref $varName"
runs[[5]] = "python calcMetricCatchments.py 4 2 24 8 $lag $model $ref $varName"
runs[[6]] = "python calcMetricCatchments.py 5 2 24 10 $lag $model $ref $varName"
runs[[7]] = "python calcMetricCatchments.py 6 2 24 12 $lag $model $ref $varName"
runs[[8]] = "python calcMetricCatchments.py 7 2 24 14 $lag $model $ref $varName"
runs[[9]] = "python calcMetricCatchments.py 8 2 24 16 $lag $model $ref $varName"
runs[[10]] = "python calcMetricCatchments.py 9 2 24 18 $lag $model $ref $varName"
runs[[11]] = "python calcMetricCatchments.py 10 2 24 20 $lag $model $ref $varName"
runs[[12]] = "python calcMetricCatchments.py 11 2 24 22 $lag $model $ref $varName"
runs[[13]] = "python calcMetricCatchments.py 12 2 24 24 $lag $model $ref $varName"

#runs[[1]] ="python calcMetrics.py 0 1 24 1 $lag $model $ref $varName"
#runs[[2]] = "python calcMetrics.py 1 2 24 2 $lag $model $ref $varName"
#runs[[3]] = "python calcMetrics.py 2 2 24 4 $lag $model $ref $varName"
#runs[[4]] = "python calcMetrics.py 3 2 24 6 $lag $model $ref $varName"
#runs[[5]] = "python calcMetrics.py 4 2 24 8 $lag $model $ref $varName"
#runs[[6]] = "python calcMetrics.py 5 2 24 10 $lag $model $ref $varName"
#runs[[7]] = "python calcMetrics.py 6 2 24 12 $lag $model $ref $varName"
#runs[[8]] = "python calcMetrics.py 7 2 24 14 $lag $model $ref $varName"
#runs[[9]] = "python calcMetrics.py 8 2 24 16 $lag $model $ref $varName"
#runs[[10]] = "python calcMetrics.py 9 2 24 18 $lag $model $ref $varName"
#runs[[11]] = "python calcMetrics.py 10 2 24 20 $lag $model $ref $varName"
#runs[[12]] = "python calcMetrics.py 11 2 24 22 $lag $model $ref $varName"
#runs[[13]] = "python calcMetrics.py 12 2 24 24 $lag $model $ref $varName"

modelS = c("WeightedEqual")#"CanCM3","CanCM4","FLOR", "EPS", "CCSM")  #c("Weighted", "WeightedEqual") #"CanCM3","CanCM4","FLOR", "EPS")
varNameS = c("evap") #,"soilMoistureLow","soilMoistureUp") #discharge") #, "groundwater","soilMoistureLow","soilMoistureUp")
refS = c("PGF")#,"CFS")
lagS = c(0:11)

jobDir = "/tigress/nwanders/Scripts/hydroSeasonal/jobs/"
master = ""

for(m in 1:length(modelS)){
  model = modelS[m]
  for(varName in varNameS){
    for(ref in refS){
      for(lag in lagS){
        jobFile = paste(jobDir,model,"_",varName,"_",ref,"_",lag,".sh",sep="")
        write.table(start, jobFile, col.names=FALSE, row.names=FALSE, quote=FALSE)
        setting=paste("model=",model,"\n","varName=",varName,"\n","ref=",ref,"\n","lag=",as.character(lag),"\n", sep="")
        write.table(setting, jobFile, col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
        write.table(runs[[1]], jobFile, col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
        #for(i in 2:(13-lag)){
        #  write.table(runs[[i]], jobFile, col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
        #}
        master = paste(master,"sbatch ", paste(jobFile,"\n",sep=""),sep="")
      }
    }
  }
}

write.table(master, paste(jobDir,"result_master.sh",sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE)

