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

file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_int.dat")

print(paste(sep="","Loading ",file))
int <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))
print(paste(sep="",file," loading"))

int <- int[order(int$Trajectory),]

int <- as_tibble(int)

n <- nchar(int$Value[1])
maxfix <- max(int$Fixation)
maxT <- max(int$Trajectory)
trnas <- unique(int$Number[which(int$Molecule == "tRNA")])
aarss <- unique(int$Number[which(int$Molecule == "aaRS")])
ntrna <- length(trnas)
naars <- length(aarss)

mask <- int[which(int$Type == "Mask"),]

mask <- mask %>% mutate(V.dec = strtoi(Value,2))

weight <- function(x){
    return(sum(as.integer(intToBits(x))))
}


mask <- mask %>% select(-Type) %>% select(-Value)

mask <- cbind(mask,sapply(mask$V.dec,weight))
colnames(mask)[dim(mask)[2]] <- "Weight"

tempmat <- matrix(mask$Weight,ncol=maxT+1)

print(paste(sep="","Finding weight quantiles."))
qmin <- apply(tempmat,1,min)
q0 <- apply(tempmat,1,quantile,probs=0.05)
q1 <- apply(tempmat,1,quantile,probs=0.25)
q2 <- apply(tempmat,1,quantile,probs=0.5)
q3 <- apply(tempmat,1,quantile,probs=0.75)
q4 <- apply(tempmat,1,quantile,probs=0.95)
qmax <- apply(tempmat,1,max)

print(paste(sep="","Making Q."))
Q <- data.frame(mask[which(mask$Trajectory == 0),2:4],qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))

colnames(Q)[c(11:dim(Q)[2])] <- c("phi","Mu","Rate","Constant")

for(Mu in MU[-1]){
    file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_int.dat")
    
    print(paste(sep="","Loading ",file))
    int <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))
    print(paste(sep="",file," loading"))
    
    int <- int[order(int$Trajectory),]
    
    int <- as_tibble(int)
    
    n <- nchar(int$Value[1])
    maxfix <- max(int$Fixation)
    maxT <- max(int$Trajectory)
    trnas <- unique(int$Number[which(int$Molecule == "tRNA")])
    aarss <- unique(int$Number[which(int$Molecule == "aaRS")])
    ntrna <- length(trnas)
    naars <- length(aarss)
    
    mask <- int[which(int$Type == "Mask"),]
    
    mask <- mask %>% mutate(V.dec = strtoi(Value,2))

    mask <- mask %>% select(-Type) %>% select(-Value)

    mask <- cbind(mask,sapply(mask$V.dec,weight))
    colnames(mask)[dim(mask)[2]] <- "Weight"
    
    tempmat <- matrix(mask$Weight,ncol=maxT+1)
    
    print(paste(sep="","Finding weight quantiles."))
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    print(paste(sep="","Making Q."))
    temp <- data.frame(mask[which(mask$Trajectory == 0),2:4],qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
    
    colnames(temp)[c(11:dim(temp)[2])] <- c("phi","Mu","Rate","Constant")

    Q <- rbind(Q,temp)
}

for(constant in CONSTANT[-1]){
    for(Mu in MU){
        file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_int.dat")
        print(paste(sep="","Loading ",file))
        int <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))
        print(paste(sep="",file," loading"))
        
        int <- int[order(int$Trajectory),]
        
        int <- as_tibble(int)
        
        n <- nchar(int$Value[1])
        maxfix <- max(int$Fixation)
        maxT <- max(int$Trajectory)
        trnas <- unique(int$Number[which(int$Molecule == "tRNA")])
        aarss <- unique(int$Number[which(int$Molecule == "aaRS")])
        ntrna <- length(trnas)
        naars <- length(aarss)
        
        mask <- int[which(int$Type == "Mask"),]
        
        mask <- mask %>% mutate(V.dec = strtoi(Value,2))
        
        mask <- mask %>% select(-Type) %>% select(-Value)
        
        mask <- cbind(mask,sapply(mask$V.dec,weight))
        colnames(mask)[dim(mask)[2]] <- "Weight"
        
        tempmat <- matrix(mask$Weight,ncol=maxT+1)
        
        print(paste(sep="","Finding weight quantiles."))
        qmin <- apply(tempmat,1,min)
        q0 <- apply(tempmat,1,quantile,probs=0.05)
        q1 <- apply(tempmat,1,quantile,probs=0.25)
        q2 <- apply(tempmat,1,quantile,probs=0.5)
        q3 <- apply(tempmat,1,quantile,probs=0.75)
        q4 <- apply(tempmat,1,quantile,probs=0.95)
        qmax <- apply(tempmat,1,max)
        
        print(paste(sep="","Making Q."))
        temp <- data.frame(mask[which(mask$Trajectory == 0),2:4],qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
        
        colnames(temp)[c(11:dim(temp)[2])] <- c("phi","Mu","Rate","Constant")
        
        Q <- rbind(Q,temp)
    }
}

rate <- RATE[2]
constant <- CONSTANT[1]

for(Mu in MU){
    file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_int.dat")
    print(paste(sep="","Loading ",file))
    int <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))
    print(paste(sep="",file," loading"))
    
    int <- int[order(int$Trajectory),]
    
    int <- as_tibble(int)
    
    n <- nchar(int$Value[1])
    maxfix <- max(int$Fixation)
    maxT <- max(int$Trajectory)
    trnas <- unique(int$Number[which(int$Molecule == "tRNA")])
    aarss <- unique(int$Number[which(int$Molecule == "aaRS")])
    ntrna <- length(trnas)
    naars <- length(aarss)
    
    mask <- int[which(int$Type == "Mask"),]
    
    mask <- mask %>% mutate(V.dec = strtoi(Value,2))

    mask <- mask %>% select(-Type) %>% select(-Value)

    mask <- cbind(mask,sapply(mask$V.dec,weight))
    colnames(mask)[dim(mask)[2]] <- "Weight"
    
    tempmat <- matrix(mask$Weight,ncol=maxT+1)
    
    print(paste(sep="","Finding weight quantiles."))
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    print(paste(sep="","Making Q."))
    temp <- data.frame(mask[which(mask$Trajectory == 0),2:4],qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
    
    colnames(temp)[c(11:dim(temp)[2])] <- c("phi","Mu","Rate","Constant")

    Q <- rbind(Q,temp)
}

for(constant in CONSTANT[-1]){
    for(Mu in MU){
        file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_int.dat")
        print(paste(sep="","Loading ",file))
        int <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))
        print(paste(sep="",file," loading"))
        
        int <- int[order(int$Trajectory),]
        
        int <- as_tibble(int)
        
        n <- nchar(int$Value[1])
        maxfix <- max(int$Fixation)
        maxT <- max(int$Trajectory)
        trnas <- unique(int$Number[which(int$Molecule == "tRNA")])
        aarss <- unique(int$Number[which(int$Molecule == "aaRS")])
        ntrna <- length(trnas)
        naars <- length(aarss)
        
        mask <- int[which(int$Type == "Mask"),]
        
        mask <- mask %>% mutate(V.dec = strtoi(Value,2))
        
        mask <- mask %>% select(-Type) %>% select(-Value)
        
        mask <- cbind(mask,sapply(mask$V.dec,weight))
        colnames(mask)[dim(mask)[2]] <- "Weight"
        
        tempmat <- matrix(mask$Weight,ncol=maxT+1)
        
        print(paste(sep="","Finding weight quantiles."))
        qmin <- apply(tempmat,1,min)
        q0 <- apply(tempmat,1,quantile,probs=0.05)
        q1 <- apply(tempmat,1,quantile,probs=0.25)
        q2 <- apply(tempmat,1,quantile,probs=0.5)
        q3 <- apply(tempmat,1,quantile,probs=0.75)
        q4 <- apply(tempmat,1,quantile,probs=0.95)
        qmax <- apply(tempmat,1,max)
        
        print(paste(sep="","Making Q."))
        temp <- data.frame(mask[which(mask$Trajectory == 0),2:4],qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
        
        colnames(temp)[c(11:dim(temp)[2])] <- c("phi","Mu","Rate","Constant")
        
        Q <- rbind(Q,temp)
    }
}

for(phi in PHI[-1]){
    for(constant in CONSTANT){
        for(rte in RATE){
            for(Mu in MU){
                file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_int.dat")
                print(paste(sep="","Loading ",file))
                int <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))
                print(paste(sep="",file," loading"))
                
                int <- int[order(int$Trajectory),]
                
                int <- as_tibble(int)
                
                n <- nchar(int$Value[1])
                maxfix <- max(int$Fixation)
                maxT <- max(int$Trajectory)
                trnas <- unique(int$Number[which(int$Molecule == "tRNA")])
                aarss <- unique(int$Number[which(int$Molecule == "aaRS")])
                ntrna <- length(trnas)
                naars <- length(aarss)
                
                mask <- int[which(int$Type == "Mask"),]
                
                mask <- mask %>% mutate(V.dec = strtoi(Value,2))
                
                mask <- mask %>% select(-Type) %>% select(-Value)
                
                mask <- cbind(mask,sapply(mask$V.dec,weight))
                colnames(mask)[dim(mask)[2]] <- "Weight"
                
                tempmat <- matrix(mask$Weight,ncol=maxT+1)
                
                print(paste(sep="","Finding weight quantiles."))
                qmin <- apply(tempmat,1,min)
                q0 <- apply(tempmat,1,quantile,probs=0.05)
                q1 <- apply(tempmat,1,quantile,probs=0.25)
                q2 <- apply(tempmat,1,quantile,probs=0.5)
                q3 <- apply(tempmat,1,quantile,probs=0.75)
                q4 <- apply(tempmat,1,quantile,probs=0.95)
                qmax <- apply(tempmat,1,max)
                
                print(paste(sep="","Making Q."))
                temp <- data.frame(mask[which(mask$Trajectory == 0),2:4],qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q2)),rep(paste(sep="","Mu = ",Mu),length(q2)),rep(gsub("_"," ",rte),length(q2)),rep(paste(sep="","constant = ",constant),length(q2)))
                
                colnames(temp)[c(11:dim(temp)[2])] <- c("phi","Mu","Rate","Constant")
                
                Q <- rbind(Q,temp)
                
            }
        }
    }
}

write.table(Q,file=paste(sep="","atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"_phi",paste(collapse="-",PHI),"_Mu",paste(collapse="-",MU),"_",codon.space,ialization,"_rate_dep_and_indep_constants",paste(collapse="-",CONSTANT),"_uniform_aas_",gsub(" ","_",proofreading),"_mask_weight.dat"),row.names=FALSE)
