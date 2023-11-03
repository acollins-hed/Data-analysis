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

file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_code.dat")

print(paste(sep="","Loading file ",file))
code <- read.table(file,header=TRUE,row.names=NULL)
print(paste(sep="",file," loaded."))

maxfix <- max(code$Fixation)
maxT <- max(code$Trajectory)

aarss <- unique(code$aaRS)
trnas <- unique(code$tRNA)
naars <- length(aarss)
ntrna <- length(trnas)
maxfix <- max(code$Fixation)
maxT <- max(code$Trajectory)

print(paste(sep="","Finding each tRNA degree."))
tdegree <- apply(matrix(code$Match,nrow=naars),2,sum)

code <- as_tibble(code)

print(paste(sep="","Pivot wider."))
temp <- code %>% select(-Prob_Interaction) %>% select(-kd) %>% pivot_wider(names_from=tRNA,values_from=Match)

print(paste(sep="","Finding each aaRS degree."))
adegree <- apply(temp[,-1:-3],1,sum)

tempmat <- matrix(adegree,ncol=maxT+1)

print(paste(sep="","Finding quantiles for aaRSs."))
qmin <- apply(tempmat,1,min)
q0 <- apply(tempmat,1,quantile,probs=0.05)
q1 <- apply(tempmat,1,quantile,probs=0.25)
q2 <- apply(tempmat,1,quantile,probs=0.5)
q3 <- apply(tempmat,1,quantile,probs=0.75)
q4 <- apply(tempmat,1,quantile,probs=0.95)
qmax <- apply(tempmat,1,max)

print(paste(sep="","Making Qa."))
Qa <- data.frame(rep(0:maxfix,each=naars),rep(aarss,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))

colnames(Qa)[c(1:2,10:dim(Qa)[2])] <- c("Fixation","aaRS","phi","Mu","Rate","Constant")

tempmat <- matrix(tdegree,ncol=maxT+1)

print(paste(sep="","Finding quantiles tRNAs."))
qmin <- apply(tempmat,1,min)
q0 <- apply(tempmat,1,quantile,probs=0.05)
q1 <- apply(tempmat,1,quantile,probs=0.25)
q2 <- apply(tempmat,1,quantile,probs=0.5)
q3 <- apply(tempmat,1,quantile,probs=0.75)
q4 <- apply(tempmat,1,quantile,probs=0.95)
qmax <- apply(tempmat,1,max)

print(paste(sep="","Making Qt."))
Qt <- data.frame(rep(0:maxfix,each=ntrna),rep(trnas,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))

colnames(Qt)[c(1:2,10:dim(Qt)[2])] <- c("Fixation","tRNA","phi","Mu","Rate","Constant")

for(Mu in MU[-1]){
    file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_code.dat")

    print(paste(sep="","Loading file ",file))
    code <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," loaded."))
    
    maxfix <- max(code$Fixation)
    maxT <- max(code$Trajectory)
    
    aarss <- unique(code$aaRS)
    trnas <- unique(code$tRNA)
    naars <- length(aarss)
    ntrna <- length(trnas)
    maxfix <- max(code$Fixation)
    maxT <- max(code$Trajectory)
    
    print(paste(sep="","Finding each tRNA degree."))
    tdegree <- apply(matrix(code$Match,nrow=naars),2,sum)
    
    code <- as_tibble(code)
    
    print(paste(sep="","Pivot wider."))
    temp <- code %>% select(-Prob_Interaction) %>% select(-kd) %>% pivot_wider(names_from=tRNA,values_from=Match)
    
    print(paste(sep="","Finding each aaRS degree."))
    adegree <- apply(temp[,-1:-3],1,sum)
    
    tempmat <- matrix(adegree,ncol=maxT+1)
    
    print(paste(sep="","Finding quantiles for aaRSs."))
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    print(paste(sep="","Making Qa."))
    atemp <- data.frame(rep(0:maxfix,each=naars),rep(aarss,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
    
    colnames(atemp)[c(1:2,10:dim(atemp)[2])] <- c("Fixation","aaRS","phi","Mu","Rate","Constant")

    Qa <- rbind(Qa,atemp)
    
    tempmat <- matrix(tdegree,ncol=maxT+1)
    
    print(paste(sep="","Finding quantiles tRNAs."))
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    print(paste(sep="","Making Qt."))
    ttemp <- data.frame(rep(0:maxfix,each=ntrna),rep(trnas,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))

    colnames(ttemp)[c(1:2,10:dim(ttemp)[2])] <- c("Fixation","tRNA","phi","Mu","Rate","Constant")

    Qt <- rbind(Qt,ttemp)
}

for(constant in CONSTANT[-1]){
    for(Mu in MU){
        file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_code.dat")

        print(paste(sep="","Loading file ",file))
        code <- read.table(file,header=TRUE,row.names=NULL)
        print(paste(sep="",file," loaded."))
        
        maxfix <- max(code$Fixation)
        maxT <- max(code$Trajectory)
        
        aarss <- unique(code$aaRS)
        trnas <- unique(code$tRNA)
        naars <- length(aarss)
        ntrna <- length(trnas)
        maxfix <- max(code$Fixation)
        maxT <- max(code$Trajectory)
        
        print(paste(sep="","Finding each tRNA degree."))
        tdegree <- apply(matrix(code$Match,nrow=naars),2,sum)
        
        code <- as_tibble(code)
        
        print(paste(sep="","Pivot wider."))
        temp <- code %>% select(-Prob_Interaction) %>% select(-kd) %>% pivot_wider(names_from=tRNA,values_from=Match)
        
        print(paste(sep="","Finding each aaRS degree."))
        adegree <- apply(temp[,-1:-3],1,sum)
        
        tempmat <- matrix(adegree,ncol=maxT+1)
        
        print(paste(sep="","Finding quantiles for aaRSs."))
        qmin <- apply(tempmat,1,min)
        q0 <- apply(tempmat,1,quantile,probs=0.05)
        q1 <- apply(tempmat,1,quantile,probs=0.25)
        q2 <- apply(tempmat,1,quantile,probs=0.5)
        q3 <- apply(tempmat,1,quantile,probs=0.75)
        q4 <- apply(tempmat,1,quantile,probs=0.95)
        qmax <- apply(tempmat,1,max)
        
        print(paste(sep="","Making Qa."))
        atemp <- data.frame(rep(0:maxfix,each=naars),rep(aarss,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
        
        colnames(atemp)[c(1:2,10:dim(atemp)[2])] <- c("Fixation","aaRS","phi","Mu","Rate","Constant")
        
        Qa <- rbind(Qa,atemp)
        
        tempmat <- matrix(tdegree,ncol=maxT+1)
        
        print(paste(sep="","Finding quantiles tRNAs."))
        qmin <- apply(tempmat,1,min)
        q0 <- apply(tempmat,1,quantile,probs=0.05)
        q1 <- apply(tempmat,1,quantile,probs=0.25)
        q2 <- apply(tempmat,1,quantile,probs=0.5)
        q3 <- apply(tempmat,1,quantile,probs=0.75)
        q4 <- apply(tempmat,1,quantile,probs=0.95)
        qmax <- apply(tempmat,1,max)
        
        print(paste(sep="","Making Qt."))
        ttemp <- data.frame(rep(0:maxfix,each=ntrna),rep(trnas,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))

        colnames(ttemp)[c(1:2,10:dim(ttemp)[2])] <- c("Fixation","tRNA","phi","Mu","Rate","Constant")
                
        Qt <- rbind(Qt,ttemp)
        
    }
}

rte <- RATE[2]
constant <- CONSTANT[1]

for(Mu in MU){
    file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_code.dat")
    
    print(paste(sep="","Loading file ",file))
    code <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," loaded."))
    
    maxfix <- max(code$Fixation)
    maxT <- max(code$Trajectory)
    
    aarss <- unique(code$aaRS)
    trnas <- unique(code$tRNA)
    naars <- length(aarss)
    ntrna <- length(trnas)
    maxfix <- max(code$Fixation)
    maxT <- max(code$Trajectory)
    
    print(paste(sep="","Finding each tRNA degree."))
    tdegree <- apply(matrix(code$Match,nrow=naars),2,sum)
    
    code <- as_tibble(code)
    
    print(paste(sep="","Pivot wider."))
    temp <- code %>% select(-Prob_Interaction) %>% select(-kd) %>% pivot_wider(names_from=tRNA,values_from=Match)
    
    print(paste(sep="","Finding each aaRS degree."))
    adegree <- apply(temp[,-1:-3],1,sum)
    
    tempmat <- matrix(adegree,ncol=maxT+1)
    
    print(paste(sep="","Finding quantiles for aaRSs."))
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    print(paste(sep="","Making Qa."))
    atemp <- data.frame(rep(0:maxfix,each=naars),rep(aarss,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
    
    colnames(atemp)[c(1:2,10:dim(atemp)[2])] <- c("Fixation","aaRS","phi","Mu","Rate","Constant")
    
    Qa <- rbind(Qa,atemp)
    
    tempmat <- matrix(tdegree,ncol=maxT+1)
    
    print(paste(sep="","Finding quantiles tRNAs."))
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    print(paste(sep="","Making Qt."))
    ttemp <- data.frame(rep(0:maxfix,each=ntrna),rep(trnas,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))

    colnames(ttemp)[c(1:2,10:dim(ttemp)[2])] <- c("Fixation","tRNA","phi","Mu","Rate","Constant")
    
    Qt <- rbind(Qt,ttemp)
}

for(constant in CONSTANT[-1]){
    for(Mu in MU){
        file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_code.dat")

        print(paste(sep="","Loading file ",file))
        code <- read.table(file,header=TRUE,row.names=NULL)
        print(paste(sep="",file," loaded."))
        
        maxfix <- max(code$Fixation)
        maxT <- max(code$Trajectory)
        
        aarss <- unique(code$aaRS)
        trnas <- unique(code$tRNA)
        naars <- length(aarss)
        ntrna <- length(trnas)
        maxfix <- max(code$Fixation)
        maxT <- max(code$Trajectory)
        
        print(paste(sep="","Finding each tRNA degree."))
        tdegree <- apply(matrix(code$Match,nrow=naars),2,sum)
        
        code <- as_tibble(code)
        
        print(paste(sep="","Pivot wider."))
        temp <- code %>% select(-Prob_Interaction) %>% select(-kd) %>% pivot_wider(names_from=tRNA,values_from=Match)
        
        print(paste(sep="","Finding each aaRS degree."))
        adegree <- apply(temp[,-1:-3],1,sum)
        
        tempmat <- matrix(adegree,ncol=maxT+1)
        
        print(paste(sep="","Finding quantiles for aaRSs."))
        qmin <- apply(tempmat,1,min)
        q0 <- apply(tempmat,1,quantile,probs=0.05)
        q1 <- apply(tempmat,1,quantile,probs=0.25)
        q2 <- apply(tempmat,1,quantile,probs=0.5)
        q3 <- apply(tempmat,1,quantile,probs=0.75)
        q4 <- apply(tempmat,1,quantile,probs=0.95)
        qmax <- apply(tempmat,1,max)
        
        print(paste(sep="","Making Qa."))
        atemp <- data.frame(rep(0:maxfix,each=naars),rep(aarss,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
        
        colnames(atemp)[c(1:2,10:dim(atemp)[2])] <- c("Fixation","aaRS","phi","Mu","Rate","Constant")
        
        Qa <- rbind(Qa,atemp)
        
        tempmat <- matrix(tdegree,ncol=maxT+1)
        
        print(paste(sep="","Finding quantiles tRNAs."))
        qmin <- apply(tempmat,1,min)
        q0 <- apply(tempmat,1,quantile,probs=0.05)
        q1 <- apply(tempmat,1,quantile,probs=0.25)
        q2 <- apply(tempmat,1,quantile,probs=0.5)
        q3 <- apply(tempmat,1,quantile,probs=0.75)
        q4 <- apply(tempmat,1,quantile,probs=0.95)
        qmax <- apply(tempmat,1,max)
        
        print(paste(sep="","Making Qt."))
        ttemp <- data.frame(rep(0:maxfix,each=ntrna),rep(trnas,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))

        colnames(ttemp)[c(1:2,10:dim(ttemp)[2])] <- c("Fixation","tRNA","phi","Mu","Rate","Constant")
        
        Qt <- rbind(Qt,ttemp)        
    }
}

for(phi in PHI[-1]){
    for(constant in CONSTANT){
        for(rte in RATE){
            for(Mu in MU){
                file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_code.dat")
                
                print(paste(sep="","Loading file ",file))
                code <- read.table(file,header=TRUE,row.names=NULL)
                print(paste(sep="",file," loaded."))
                
                maxfix <- max(code$Fixation)
                maxT <- max(code$Trajectory)
                
                aarss <- unique(code$aaRS)
                trnas <- unique(code$tRNA)
                naars <- length(aarss)
                ntrna <- length(trnas)
                maxfix <- max(code$Fixation)
                maxT <- max(code$Trajectory)
                
                print(paste(sep="","Finding each tRNA degree."))
                tdegree <- apply(matrix(code$Match,nrow=naars),2,sum)
                
                code <- as_tibble(code)
                
                print(paste(sep="","Pivot wider."))
                temp <- code %>% select(-Prob_Interaction) %>% select(-kd) %>% pivot_wider(names_from=tRNA,values_from=Match)
                
                print(paste(sep="","Finding each aaRS degree."))
                adegree <- apply(temp[,-1:-3],1,sum)
                
                tempmat <- matrix(adegree,ncol=maxT+1)
                
                print(paste(sep="","Finding quantiles for aaRSs."))
                qmin <- apply(tempmat,1,min)
                q0 <- apply(tempmat,1,quantile,probs=0.05)
                q1 <- apply(tempmat,1,quantile,probs=0.25)
                q2 <- apply(tempmat,1,quantile,probs=0.5)
                q3 <- apply(tempmat,1,quantile,probs=0.75)
                q4 <- apply(tempmat,1,quantile,probs=0.95)
                qmax <- apply(tempmat,1,max)
                
                print(paste(sep="","Making Qa."))
                atemp <- data.frame(rep(0:maxfix,each=naars),rep(aarss,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
                
                colnames(atemp)[c(1:2,10:dim(atemp)[2])] <- c("Fixation","aaRS","phi","Mu","Rate","Constant")
                
                Qa <- rbind(Qa,atemp)
                
                tempmat <- matrix(tdegree,ncol=maxT+1)
                
                print(paste(sep="","Finding quantiles tRNAs."))
                qmin <- apply(tempmat,1,min)
                q0 <- apply(tempmat,1,quantile,probs=0.05)
                q1 <- apply(tempmat,1,quantile,probs=0.25)
                q2 <- apply(tempmat,1,quantile,probs=0.5)
                q3 <- apply(tempmat,1,quantile,probs=0.75)
                q4 <- apply(tempmat,1,quantile,probs=0.95)
                qmax <- apply(tempmat,1,max)
                
                print(paste(sep="","Making Qt."))
                ttemp <- data.frame(rep(0:maxfix,each=ntrna),rep(trnas,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))

                colnames(ttemp)[c(1:2,10:dim(ttemp)[2])] <- c("Fixation","tRNA","phi","Mu","Rate","Constant")
        
                Qt <- rbind(Qt,ttemp)
            }
        }
    }
}

write.table(Qa,file=paste(sep="","atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"_phi",paste(collapse="-",PHI),"_Mu",paste(collapse="-",MU),"_",codon.space,ialization,"_rate_dep_and_indep_constants",paste(collapse="-",CONSTANT),"_uniform_aas_",gsub(" ","_",proofreading),"_aaRS_degree.dat"),row.names=FALSE)

write.table(Qt,file=paste(sep="","atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"_phi",paste(collapse="-",PHI),"_Mu",paste(collapse="-",MU),"_",codon.space,ialization,"_rate_dep_and_indep_constants",paste(collapse="-",CONSTANT),"_uniform_aas_",gsub(" ","_",proofreading),"_tRNA_degree.dat"),row.names=FALSE)
