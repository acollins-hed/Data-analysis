require(tidyverse)
codon.frequencies <- c("all","frequencies are 10")
rate <- "rate dep"
initialization <- "all_0s"
init <- "all 0s"
codon.space <- "Cring"
cs <- "codon ring"
phi <- 0.5
k <- 4
n <- 8
naars <- 5
ntrna <- 5

MU <- c("0.0001","0.001","0.005","0.01","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2")

Mu <- MU[1]

file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_prob.dat")

print(paste(sep="","Loading ",file))
prob <- read.table(file,header=TRUE,row.names=NULL)
print(paste(sep="",file," Loaded"))

AAs <- unique(prob$Amino_Acid)
aas <- length(AAs)
Sts <- unique(prob$Site_Type)
sts <- length(Sts)
p.of.s <- 1/sts

maxfix <- max(prob$Fixation)
maxT <- max(prob$Trajectory)

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

tempmat <- matrix(p.of.a$EX.term,ncol=naars)

EX <- apply(tempmat,1,sum)

tempmat <- matrix(p.of.a$EX.sq.term,ncol=naars)

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

Q <- as.data.frame(cbind(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","Mu = ",Mu),length(q0))))

colnames(Q)[c(1,length(colnames(Q)))] <- c("Fixation","Mu")

for(Mu in MU[-1]){

    file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_prob.dat")

    print(paste(sep="","Loading ",file))
    prob <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," Loaded"))
    
    AAs <- unique(prob$Amino_Acid)
    aas <- length(AAs)
    Sts <- unique(prob$Site_Type)
    sts <- length(Sts)
    p.of.s <- 1/sts
    
    maxfix <- max(prob$Fixation)
    maxT <- max(prob$Trajectory)
    
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
    
    tempmat <- matrix(p.of.a$EX.term,ncol=naars)
    
    EX <- apply(tempmat,1,sum)
    
    tempmat <- matrix(p.of.a$EX.sq.term,ncol=naars)
    
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
    
    temp <- as.data.frame(cbind(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","Mu = ",Mu),length(q0))))
    
    colnames(temp)[c(1,length(colnames(temp)))] <- c("Fixation","Mu")

    Q <- rbind(Q,temp)
}

Q$Fixation <- as.integer(Q$Fixation)
Q$qmin <- as.numeric(Q$qmin)
Q$q0 <- as.numeric(Q$q0)
Q$q1 <- as.numeric(Q$q1)
Q$q2 <- as.numeric(Q$q2)
Q$q3 <- as.numeric(Q$q3)
Q$q4 <- as.numeric(Q$q4)
Q$qmax <- as.numeric(Q$qmax)

image.name <- paste(sep="","P_Of_AA_St.Dev._Vs_Mu_S",sts,"AA",aas,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- bquote(2*sigma~"[v("~alpha~")] "~.(paste(sep="","Over ",maxT+1))~"Trajectories,"~.(paste(sep="",maxfix+1," Fixations Each")))
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",naars,", T = ",ntrna,", S = ",sts,", n = ",n,", k = ",k,",")
image.subtitles[2] <- paste(sep="","= ",phi,", ",cs,", ",rate,", Initialized ",init,",")

png(image.name,height=1350,width=2200)
print(ggplot(Q,aes(x=Fixation,y=q2))+geom_line(aes(color="black"))+geom_ribbon(data=Q,aes(ymin=q1,ymax=q3),color="green",alpha=0.3)+geom_ribbon(data=Q,aes(ymin=q0,ymax=q4),color="purple",alpha=0.2)+geom_ribbon(data=Q,aes(ymin=qmin,ymax=qmax),color="red",alpha=0.1)+facet_grid(cols=vars(Mu))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+theme(legend.position="bottom",text=element_text(size=30),strip.text=element_text(size=25),axis.text.x=element_text(size=15))+ylab(bquote(2*sigma~"[v("~alpha~")]"))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()

write.table(Q,"A5T5S5_phi0.5_all0_rate_dep_cring_P_Of_AA_St.Dev_Quantiles.dat")
