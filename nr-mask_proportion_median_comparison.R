require(tidyverse)

codon.frequencies <- c("all","frequencies are 10")
RATE <- c("rate_dep","rate_indep")
rate <- gsub("_"," ",RATE[1])
initialization <- "all_1s"
init <- gsub("_"," ",initialization)
ialization <- gsub("s","",gsub("_","",initialization))
codon.space <- "Cring"
cs <- "codon ring"
s <- 5
kmin <- 2
A <- 5
T <- 5
n <- 8
proofreading <- "proofreading"

#PHI <-c("0.01","0.1","0.25","0.5","0.75","0.9","0.99")
#PHI <-c("0.5")
PHI <-c("0.01","0.1","0.25","0.5","0.75","0.9")
#MU <- c("0.0001","0.001","0.01","0.1","0.2")
#MU <- c("0.01","0.1","0.2")
MU <- c("0.0001","0.001","0.01","0.1","0.2","0.3")
#CONSTANT <- c("10","100","1000","10000","100000")
#CONSTANT <- c("1e08","1e09","1e10","1e11","1e12")
CONSTANT <- c("1e-09")

phi <- PHI[1]
Mu <- MU[1]
rte <- RATE[1]
constant <- CONSTANT[1]

#directory <- "/media/sciship/7D6B-3CD0/Data/a4t4n8k4C2/"
directory <- "/media/sciship/Seagate\ Expansion\ Drive/New_Rate_Function/a5t5n8kmin2Cring_proofreading/"

file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_traj.dat")

print(paste(sep="","Loading file ",file))
traj <- read.table(file,header=TRUE,row.names=NULL)
print(paste(sep="",file," loaded."))

maxfix <- max(traj$Fixation)
maxT <- max(traj$Trajectory)

tempmat <- matrix(traj$Percent_On,nrow=maxfix+1)

Median <- apply(tempmat,1,median)

Q <- data.frame(0:maxfix,Median,rep(paste(sep="","Mu = ",Mu),length(Median)),rep(paste(sep="","phi = ",phi),length(Median)),rep(gsub("_"," ",rte),length(Median)),rep(paste(sep="","constant = ",constant),length(Median)))
colnames(Q) <- c("Fixation","Median","Mu","phi","Rate","Constant")

for(Mu in MU[-1]){
    file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_traj.dat")

    print(paste(sep="","Loading file ",file))
    traj <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," loaded."))
    
    maxfix <- max(traj$Fixation)
    maxT <- max(traj$Trajectory)
    
    tempmat <- matrix(traj$Percent_On,nrow=maxfix+1)
    
    Median <- apply(tempmat,1,median)

    temp <- data.frame(0:maxfix,Median,rep(paste(sep="","Mu = ",Mu),length(Median)),rep(paste(sep="","phi = ",phi),length(Median)),rep(gsub("_"," ",rte),length(Median)),rep(paste(sep="","constant = ",constant),length(Median)))
    colnames(temp) <- c("Fixation","Median","Mu","phi","Rate","Constant")

    Q <- rbind(Q,temp)
}

for(constant in CONSTANT[-1]){
    for(Mu in MU){
        file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_traj.dat")
        
        print(paste(sep="","Loading file ",file))
        traj <- read.table(file,header=TRUE,row.names=NULL)
        print(paste(sep="",file," loaded."))
        
        maxfix <- max(traj$Fixation)
        maxT <- max(traj$Trajectory)
        
        tempmat <- matrix(traj$Percent_On,nrow=maxfix+1)
        
        Median <- apply(tempmat,1,median)
        
        temp <- data.frame(0:maxfix,Median,rep(paste(sep="","Mu = ",Mu),length(Median)),rep(paste(sep="","phi = ",phi),length(Median)),rep(gsub("_"," ",rte),length(Median)),rep(paste(sep="","constant = ",constant),length(Median)))
        colnames(temp) <- c("Fixation","Median","Mu","phi","Rate","Constant")
        
        Q <- rbind(Q,temp)
    }
}

rte <- RATE[2]
constant <- CONSTANT[1]

for(Mu in MU){
    file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_traj.dat")

    print(paste(sep="","Loading file ",file))
    traj <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," loaded."))
    
    maxfix <- max(traj$Fixation)
    maxT <- max(traj$Trajectory)
    
    tempmat <- matrix(traj$Percent_On,nrow=maxfix+1)
    
    Median <- apply(tempmat,1,median)

    temp <- data.frame(0:maxfix,Median,rep(paste(sep="","Mu = ",Mu),length(Median)),rep(paste(sep="","phi = ",phi),length(Median)),rep(gsub("_"," ",rte),length(Median)),rep(paste(sep="","constant = ",constant),length(Median)))
    colnames(temp) <- c("Fixation","Median","Mu","phi","Rate","Constant")

    Q <- rbind(Q,temp)
}

for(constant in CONSTANT[-1]){
    for(Mu in MU){
        file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_traj.dat")
        
        print(paste(sep="","Loading file ",file))
        traj <- read.table(file,header=TRUE,row.names=NULL)
        print(paste(sep="",file," loaded."))
        
        maxfix <- max(traj$Fixation)
        maxT <- max(traj$Trajectory)
        
        tempmat <- matrix(traj$Percent_On,nrow=maxfix+1)
        
        Median <- apply(tempmat,1,median)
        
        temp <- data.frame(0:maxfix,Median,rep(paste(sep="","Mu = ",Mu),length(Median)),rep(paste(sep="","phi = ",phi),length(Median)),rep(gsub("_"," ",rte),length(Median)),rep(paste(sep="","constant = ",constant),length(Median)))
        colnames(temp) <- c("Fixation","Median","Mu","phi","Rate","Constant")
        
        Q <- rbind(Q,temp)
    }
}


for(phi in PHI[-1]){
    for(constant in CONSTANT){
        for(rte in RATE){
            for(Mu in MU){
                file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_traj.dat")
                
                print(paste(sep="","Loading file ",file))
                traj <- read.table(file,header=TRUE,row.names=NULL)
                print(paste(sep="",file," loaded."))
                
                maxfix <- max(traj$Fixation)
                maxT <- max(traj$Trajectory)
                
                tempmat <- matrix(traj$Percent_On,nrow=maxfix+1)
                
                Median <- apply(tempmat,1,median)
                
                temp <- data.frame(0:maxfix,Median,rep(paste(sep="","Mu = ",Mu),length(Median)),rep(paste(sep="","phi = ",phi),length(Median)),rep(gsub("_"," ",rte),length(Median)),rep(paste(sep="","constant = ",constant),length(Median)))
                colnames(temp) <- c("Fixation","Median","Mu","phi","Rate","Constant")
                
                Q <- rbind(Q,temp)
            }       
        }
    }
}

write.table(Q,file=paste(sep="","atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"_phi",paste(collapse="-",PHI),"_Mu",paste(collapse="-",MU),"_",codon.space,ialization,"_rate_dep_and_indep_constants",paste(collapse="-",CONSTANT),"_uniform_aas_",gsub(" ","_",proofreading),"_mask_proportion_median_comparison.dat"),row.names=FALSE)

image.name <- paste(sep="","Proportion_ON_Median_Comparison_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"_phi",paste(collapse="-",PHI),"_Mu",paste(collapse="-",MU),"_",codon.space,"_",initialization,"_rate_dep_and_indep_constant",paste(collapse="-",CONSTANT),"_",gsub(" ","_",proofreading),".png")
image.title <- paste(sep="","Median Proportion ON Bits Over ",maxfix+1," Fixations Over ",maxT+1," Trajectories")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",A,", T = ",T,", S = ",s,", n = ",n,", kmin = ",kmin,",")
image.subtitles[2] <- paste(sep=""," ",cs,", Initialized ",init,", ",proofreading,",")

png(image.name,width=2200,height=1350)
print(ggplot(Q,aes(x=Fixation,y=Median,color=Rate))+geom_line()+facet_grid(cols=vars(Mu),rows=vars(Constant))+theme(text = element_text(size=30),strip.text.y=element_text(size=25),axis.text.y=element_text(size=20),legend.position="bottom")+ylab(bquote("Median Proportion ON Bits"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~"= 0.5,"~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2]))))
dev.off()

image.name <- paste(sep="","Proportion_ON_Median_Comparison_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"_phi",paste(collapse="-",PHI),"_Mu",paste(collapse="-",MU),"_",codon.space,"_",initialization,"_rate_dep_and_indep_constant",paste(collapse="-",CONSTANT),"_",gsub(" ","_",proofreading),".png")
image.title <- paste(sep="","Median Proportion ON Bits Over ",maxfix+1," Fixations Over ",maxT+1," Trajectories")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",A,", T = ",T,", S = ",s,", n = ",n,", kmin = ",kmin,",")
image.subtitles[2] <- paste(sep=""," ",cs,", Initialized ",init,", ",proofreading,",")

png(image.name,width=2200,height=1350)
print(ggplot(Q,aes(x=Fixation,y=Median,color=Rate))+geom_line()+facet_grid(cols=vars(Mu),rows=vars(phi))+theme(text = element_text(size=30),strip.text.y=element_text(size=25),axis.text.y=element_text(size=20),legend.position="bottom")+ylab(bquote("Median Proportion ON Bits"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~" f = 1e-09,"~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2]))))
dev.off()
