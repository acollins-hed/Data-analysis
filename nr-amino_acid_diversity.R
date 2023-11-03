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

file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_prob.dat")

print(paste(sep="","Loading file ",file))
prob <- read.table(file,header=TRUE,row.names=NULL)
print(paste(sep="",file," loaded."))

maxfix <- max(prob$Fixation)
maxT <- max(prob$Trajectory)

AAs <- unique(prob$Amino_Acid)
aas <- length(AAs)
Sts <- unique(prob$Site_Type)
sts <- length(Sts)
p.of.s <- 1/sts

p.g.s <- prob$Probability*p.of.s
p.of.a <- c(0)

print(paste(sep="","Beginning calculations for p.o.a"))
for(a in 1:aas){
    p.o.a <- apply(matrix(p.g.s[which(prob$Amino_Acid == AAs[a])],nrow=sts),2,sum)
    p.o.a <- cbind(rep(0:maxT,each=maxfix+1),rep(0:maxfix,maxT+1),rep(AAs[a],(maxfix+1)*(maxT+1)),p.o.a)
    if(a == 1){
        p.of.a <- p.o.a
    }else{
        p.of.a <- rbind(p.of.a,p.o.a)
    }
}

p.of.a <- as.data.frame(p.of.a)
colnames(p.of.a) <- c("Trajectory","Fixation","Amino_Acid","Probability")

p.of.a <- as_tibble(p.of.a)

print("Calculating First Standard Deviation Range")
p.of.a <- p.of.a %>% mutate(EX.term = Amino_Acid*Probability,EX.sq.term=Amino_Acid^2*Probability)

tempmat <- matrix(p.of.a$EX.term,ncol=A)

EX <- apply(tempmat,1,sum)

tempmat <- matrix(p.of.a$EX.sq.term,ncol=A)

EX.sq <- apply(tempmat,1,sum)

St.Dev.Range <- 2*(EX.sq - EX^2)^(0.5)

tempmat <- matrix(St.Dev.Range,nrow=maxfix+1)

print("Calculating quantiles")
qmin <- apply(tempmat,1,min)
q0 <- apply(tempmat,1,quantile,probs=0.05)
q1 <- apply(tempmat,1,quantile,probs=0.25)
q2 <- apply(tempmat,1,quantile,probs=0.5)
q3 <- apply(tempmat,1,quantile,probs=0.75)
q4 <- apply(tempmat,1,quantile,probs=0.95)
qmax <- apply(tempmat,1,max)

Q <- data.frame(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q0)),rep(paste(sep="","Mu = ",Mu),length(q0)),rep(paste(sep="","constant = ",constant),length(q0)),rep(gsub("_"," ",rte),length(q0)))

colnames(Q)[c(1,9:length(colnames(Q)))] <- c("Fixation","phi","Mu","Constant","Rate")

for(Mu in MU[-1]){
    file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_prob.dat")

    print(paste(sep="","Loading file ",file))
    prob <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," loaded."))
    
    maxfix <- max(prob$Fixation)
    maxT <- max(prob$Trajectory)
    
    AAs <- unique(prob$Amino_Acid)
    aas <- length(AAs)
    Sts <- unique(prob$Site_Type)
    sts <- length(Sts)
    p.of.s <- 1/sts
    
    p.g.s <- prob$Probability*p.of.s
    p.of.a <- c(0)
    
    print(paste(sep="","Beginning calculations for p.o.a"))
    for(a in 1:aas){
        p.o.a <- apply(matrix(p.g.s[which(prob$Amino_Acid == AAs[a])],nrow=sts),2,sum)
        p.o.a <- cbind(rep(0:maxT,each=maxfix+1),rep(0:maxfix,maxT+1),rep(AAs[a],(maxfix+1)*(maxT+1)),p.o.a)
        if(a == 1){
            p.of.a <- p.o.a
        }else{
            p.of.a <- rbind(p.of.a,p.o.a)
        }
    }
    
    p.of.a <- as.data.frame(p.of.a)
    colnames(p.of.a) <- c("Trajectory","Fixation","Amino_Acid","Probability")
    
    p.of.a <- as_tibble(p.of.a)
    
    print("Calculating First Standard Deviation Range")
    p.of.a <- p.of.a %>% mutate(EX.term = Amino_Acid*Probability,EX.sq.term=Amino_Acid^2*Probability)
    
    tempmat <- matrix(p.of.a$EX.term,ncol=A)
    
    EX <- apply(tempmat,1,sum)
    
    tempmat <- matrix(p.of.a$EX.sq.term,ncol=A)
    
    EX.sq <- apply(tempmat,1,sum)
    
    St.Dev.Range <- 2*(EX.sq - EX^2)^(0.5)
    
    tempmat <- matrix(St.Dev.Range,nrow=maxfix+1)
    
    print("Calculating quantiles")
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    temp <- data.frame(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q0)),rep(paste(sep="","Mu = ",Mu),length(q0)),rep(paste(sep="","constant = ",constant),length(q0)),rep(gsub("_"," ",rte),length(q0)))
    
    colnames(temp)[c(1,9:length(colnames(temp)))] <- c("Fixation","phi","Mu","Constant","Rate")

    print("Binding Q")
    Q <- rbind(Q,temp)
}

for(constant in CONSTANT[-1]){
    for(Mu in MU){
        file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_prob.dat")
        
        print(paste(sep="","Loading file ",file))
        prob <- read.table(file,header=TRUE,row.names=NULL)
        print(paste(sep="",file," loaded."))
        
        maxfix <- max(prob$Fixation)
        maxT <- max(prob$Trajectory)
        
        AAs <- unique(prob$Amino_Acid)
        aas <- length(AAs)
        Sts <- unique(prob$Site_Type)
        sts <- length(Sts)
        p.of.s <- 1/sts
        
        p.g.s <- prob$Probability*p.of.s
        p.of.a <- c(0)
        
        print(paste(sep="","Beginning calculations for p.o.a"))
        for(a in 1:aas){
            p.o.a <- apply(matrix(p.g.s[which(prob$Amino_Acid == AAs[a])],nrow=sts),2,sum)
            p.o.a <- cbind(rep(0:maxT,each=maxfix+1),rep(0:maxfix,maxT+1),rep(AAs[a],(maxfix+1)*(maxT+1)),p.o.a)
            if(a == 1){
                p.of.a <- p.o.a
            }else{
                p.of.a <- rbind(p.of.a,p.o.a)
            }
        }
        
        p.of.a <- as.data.frame(p.of.a)
        colnames(p.of.a) <- c("Trajectory","Fixation","Amino_Acid","Probability")
        
        p.of.a <- as_tibble(p.of.a)
        
        print("Calculating First Standard Deviation Range")
        p.of.a <- p.of.a %>% mutate(EX.term = Amino_Acid*Probability,EX.sq.term=Amino_Acid^2*Probability)
        
        tempmat <- matrix(p.of.a$EX.term,ncol=A)
        
        EX <- apply(tempmat,1,sum)
        
        tempmat <- matrix(p.of.a$EX.sq.term,ncol=A)
        
        EX.sq <- apply(tempmat,1,sum)
        
        St.Dev.Range <- 2*(EX.sq - EX^2)^(0.5)
        
        tempmat <- matrix(St.Dev.Range,nrow=maxfix+1)
        
        print("Calculating quantiles")
        qmin <- apply(tempmat,1,min)
        q0 <- apply(tempmat,1,quantile,probs=0.05)
        q1 <- apply(tempmat,1,quantile,probs=0.25)
        q2 <- apply(tempmat,1,quantile,probs=0.5)
        q3 <- apply(tempmat,1,quantile,probs=0.75)
        q4 <- apply(tempmat,1,quantile,probs=0.95)
        qmax <- apply(tempmat,1,max)
        
        temp <- data.frame(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q0)),rep(paste(sep="","Mu = ",Mu),length(q0)),rep(paste(sep="","constant = ",constant),length(q0)),rep(gsub("_"," ",rte),length(q0)))
        
        colnames(temp)[c(1,9:length(colnames(temp)))] <- c("Fixation","phi","Mu","Constant","Rate")
        
        print("Binding Q")
        Q <- rbind(Q,temp)
    }
}

rte <- RATE[2]
constant <- CONSTANT[1]

for(Mu in MU){
    file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_prob.dat")

    print(paste(sep="","Loading file ",file))
    prob <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," loaded."))
    
    maxfix <- max(prob$Fixation)
    maxT <- max(prob$Trajectory)
    
    AAs <- unique(prob$Amino_Acid)
    aas <- length(AAs)
    Sts <- unique(prob$Site_Type)
    sts <- length(Sts)
    p.of.s <- 1/sts
    
    p.g.s <- prob$Probability*p.of.s
    p.of.a <- c(0)
    
    print(paste(sep="","Beginning calculations for p.o.a"))
    for(a in 1:aas){
        p.o.a <- apply(matrix(p.g.s[which(prob$Amino_Acid == AAs[a])],nrow=sts),2,sum)
        p.o.a <- cbind(rep(0:maxT,each=maxfix+1),rep(0:maxfix,maxT+1),rep(AAs[a],(maxfix+1)*(maxT+1)),p.o.a)
        if(a == 1){
            p.of.a <- p.o.a
        }else{
            p.of.a <- rbind(p.of.a,p.o.a)
        }
    }
    
    p.of.a <- as.data.frame(p.of.a)
    colnames(p.of.a) <- c("Trajectory","Fixation","Amino_Acid","Probability")
    
    p.of.a <- as_tibble(p.of.a)
    
    print("Calculating First Standard Deviation Range")
    p.of.a <- p.of.a %>% mutate(EX.term = Amino_Acid*Probability,EX.sq.term=Amino_Acid^2*Probability)
    
    tempmat <- matrix(p.of.a$EX.term,ncol=A)
    
    EX <- apply(tempmat,1,sum)
    
    tempmat <- matrix(p.of.a$EX.sq.term,ncol=A)
    
    EX.sq <- apply(tempmat,1,sum)
    
    St.Dev.Range <- 2*(EX.sq - EX^2)^(0.5)
    
    tempmat <- matrix(St.Dev.Range,nrow=maxfix+1)
    
    print("Calculating quantiles")
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    temp <- data.frame(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q0)),rep(paste(sep="","Mu = ",Mu),length(q0)),rep(paste(sep="","constant = ",constant),length(q0)),rep(gsub("_"," ",rte),length(q0)))
    
    colnames(temp)[c(1,9:length(colnames(temp)))] <- c("Fixation","phi","Mu","Constant","Rate")

    print("Binding Q")
    Q <- rbind(Q,temp)
}

for(constant in CONSTANT[-1]){
    for(Mu in MU){
        file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_prob.dat")
        
        print(paste(sep="","Loading file ",file))
        prob <- read.table(file,header=TRUE,row.names=NULL)
        print(paste(sep="",file," loaded."))
        
        maxfix <- max(prob$Fixation)
        maxT <- max(prob$Trajectory)
        
        AAs <- unique(prob$Amino_Acid)
        aas <- length(AAs)
        Sts <- unique(prob$Site_Type)
        sts <- length(Sts)
        p.of.s <- 1/sts
        
        p.g.s <- prob$Probability*p.of.s
        p.of.a <- c(0)
        
        print(paste(sep="","Beginning calculations for p.o.a"))
        for(a in 1:aas){
            p.o.a <- apply(matrix(p.g.s[which(prob$Amino_Acid == AAs[a])],nrow=sts),2,sum)
            p.o.a <- cbind(rep(0:maxT,each=maxfix+1),rep(0:maxfix,maxT+1),rep(AAs[a],(maxfix+1)*(maxT+1)),p.o.a)
            if(a == 1){
                p.of.a <- p.o.a
            }else{
                p.of.a <- rbind(p.of.a,p.o.a)
            }
        }
        
        p.of.a <- as.data.frame(p.of.a)
        colnames(p.of.a) <- c("Trajectory","Fixation","Amino_Acid","Probability")
        
        p.of.a <- as_tibble(p.of.a)
        
        print("Calculating First Standard Deviation Range")
        p.of.a <- p.of.a %>% mutate(EX.term = Amino_Acid*Probability,EX.sq.term=Amino_Acid^2*Probability)
        
        tempmat <- matrix(p.of.a$EX.term,ncol=A)
        
        EX <- apply(tempmat,1,sum)
        
        tempmat <- matrix(p.of.a$EX.sq.term,ncol=A)
        
        EX.sq <- apply(tempmat,1,sum)
        
        St.Dev.Range <- 2*(EX.sq - EX^2)^(0.5)
        
        tempmat <- matrix(St.Dev.Range,nrow=maxfix+1)
        
        print("Calculating quantiles")
        qmin <- apply(tempmat,1,min)
        q0 <- apply(tempmat,1,quantile,probs=0.05)
        q1 <- apply(tempmat,1,quantile,probs=0.25)
        q2 <- apply(tempmat,1,quantile,probs=0.5)
        q3 <- apply(tempmat,1,quantile,probs=0.75)
        q4 <- apply(tempmat,1,quantile,probs=0.95)
        qmax <- apply(tempmat,1,max)
        
        temp <- data.frame(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q0)),rep(paste(sep="","Mu = ",Mu),length(q0)),rep(paste(sep="","constant = ",constant),length(q0)),rep(gsub("_"," ",rte),length(q0)))
        
        colnames(temp)[c(1,9:length(colnames(temp)))] <- c("Fixation","phi","Mu","Constant","Rate")
        
        print("Binding Q")
        Q <- rbind(Q,temp)
    }
}

for(phi in PHI[-1]){
    for(constant in CONSTANT){
        for(rte in RATE){
            for(Mu in MU){
                file <- paste(sep="",directory,"atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"phi",phi,"Mu",Mu,"fix1000",codon.space,ialization,"_",rte,"_",gsub(" ","_",proofreading),"_rc",constant,"_unif_aas_prob.dat")
                
                print(paste(sep="","Loading file ",file))
                prob <- read.table(file,header=TRUE,row.names=NULL)
                print(paste(sep="",file," loaded."))
                
                maxfix <- max(prob$Fixation)
                maxT <- max(prob$Trajectory)
                
                AAs <- unique(prob$Amino_Acid)
                aas <- length(AAs)
                Sts <- unique(prob$Site_Type)
                sts <- length(Sts)
                p.of.s <- 1/sts
                
                p.g.s <- prob$Probability*p.of.s
                p.of.a <- c(0)
                
                print(paste(sep="","Beginning calculations for p.o.a"))
                for(a in 1:aas){
                    p.o.a <- apply(matrix(p.g.s[which(prob$Amino_Acid == AAs[a])],nrow=sts),2,sum)
                    p.o.a <- cbind(rep(0:maxT,each=maxfix+1),rep(0:maxfix,maxT+1),rep(AAs[a],(maxfix+1)*(maxT+1)),p.o.a)
                    if(a == 1){
                        p.of.a <- p.o.a
                    }else{
                        p.of.a <- rbind(p.of.a,p.o.a)
                    }
                }
                
                p.of.a <- as.data.frame(p.of.a)
                colnames(p.of.a) <- c("Trajectory","Fixation","Amino_Acid","Probability")
                
                p.of.a <- as_tibble(p.of.a)
                
                print("Calculating First Standard Deviation Range")
                p.of.a <- p.of.a %>% mutate(EX.term = Amino_Acid*Probability,EX.sq.term=Amino_Acid^2*Probability)
                
                tempmat <- matrix(p.of.a$EX.term,ncol=A)
                
                EX <- apply(tempmat,1,sum)
                
                tempmat <- matrix(p.of.a$EX.sq.term,ncol=A)
                
                EX.sq <- apply(tempmat,1,sum)
                
                St.Dev.Range <- 2*(EX.sq - EX^2)^(0.5)
                
                tempmat <- matrix(St.Dev.Range,nrow=maxfix+1)
                
                print("Calculating quantiles")
                qmin <- apply(tempmat,1,min)
                q0 <- apply(tempmat,1,quantile,probs=0.05)
                q1 <- apply(tempmat,1,quantile,probs=0.25)
                q2 <- apply(tempmat,1,quantile,probs=0.5)
                q3 <- apply(tempmat,1,quantile,probs=0.75)
                q4 <- apply(tempmat,1,quantile,probs=0.95)
                qmax <- apply(tempmat,1,max)
                
                temp <- data.frame(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","phi = ",phi),length(q0)),rep(paste(sep="","Mu = ",Mu),length(q0)),rep(paste(sep="","constant = ",constant),length(q0)),rep(gsub("_"," ",rte),length(q0)))
                
                colnames(temp)[c(1,9:length(colnames(temp)))] <- c("Fixation","phi","Mu","Constant","Rate")
                
                print("Binding Q")
                Q <- rbind(Q,temp)
            }
        }
    }
}

write.table(Q,file=paste(sep="","atINFLATE_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"_phi",paste(collapse="-",PHI),"_Mu",paste(collapse="-",MU),"_",codon.space,ialization,"_rate_dep_and_indep_constants",paste(collapse="-",CONSTANT),"_uniform_aas_",gsub(" ","_",proofreading),"_amino_acid_diversity.dat"),row.names=FALSE)

image.name <- paste(sep="","Amino_Acid_Diversity_A",A,"T",T,"S",s,"n",n,"kmin",kmin,"_phi",paste(collapse="-",PHI),"_Mu",paste(collapse="-",MU),"_",codon.space,"_",initialization,"_rate_dep_and_indep_constant",paste(collapse="-",CONSTANT),"_",gsub(" ","_",proofreading),".png")
image.title <- paste(sep="","Amino Acid Diversity Over ",maxfix+1," Fixations Over ",maxT+1," Trajectories")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",A,", T = ",T,", S = ",s,", n = ",n,", kmin = ",kmin,",")
image.subtitles[2] <- paste(sep=""," ",cs,", Initialized ",init,", ",proofreading,",")

print(ggplot(Q,aes(x=Fixation,y=q2,color=Rate))+geom_line()+facet_grid(cols=vars(Mu),rows=vars(phi))+theme(text = element_text(size=30),strip.text.y=element_text(size=25),axis.text.y=element_text(size=20),legend.position="bottom")+ylab(bquote("Amino Acid Diversity"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~" f = 1e-09,"~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2]))))

png(image.name,height=1350,width=2200)
print(ggplot(Q,aes(x=Fixation,y=q2,color=Rate))+geom_line()+facet_grid(cols=vars(Mu),rows=vars(phi))+theme(text = element_text(size=30),strip.text.y=element_text(size=25),axis.text.y=element_text(size=20),legend.position="bottom")+ylab(bquote("Amino Acid Diversity"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~" f = 1e-09,"~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2]))))
dev.off()
