library(tidyverse)

codon.frequencies <- c("all","frequencies are 10")
rate <- "rate dep"
initialization <- "all_0s"
init <- "all 0s"
codon.space <- "Cring"
cs <- "codon ring"
s <- 5
phi <- 0.5
k <- 4
A <- 5
T <- 5
n <- 8


MU <- c("0.0001","0.001","0.005","0.01","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2")

Mu <- MU[1]

file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_traj.dat")

traj <- read.table(file,header=TRUE,row.names=NULL)

maxfix <- max(traj$Fixation)
maxT <- max(traj$Trajectory)

tempmat <- matrix(traj$Percent_On,nrow=maxfix+1)

qmin <- apply(tempmat,1,min)
q0 <- apply(tempmat,1,quantile,probs=0.05)
q1 <- apply(tempmat,1,quantile,probs=0.25)
q2 <- apply(tempmat,1,quantile,probs=0.5)
q3 <- apply(tempmat,1,quantile,probs=0.75)
q4 <- apply(tempmat,1,quantile,probs=0.95)
qmax <- apply(tempmat,1,max)


pon <- cbind(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","Mu = ",Mu),))
colnames(pon)[c(1,9)] <- c("Fixation","Mu")
pon <- as.data.frame(pon)

for(Mu in MU[-1]){

    file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_traj.dat")

    traj <- read.table(file,header=TRUE,row.names=NULL)

    maxfix <- max(traj$Fixation)
    maxT <- max(traj$Trajectory)

    tempmat <- matrix(traj$Percent_On,nrow=maxfix+1)

    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)


    temp <- cbind(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","Mu = ",Mu),))
    colnames(temp)[c(1,9)] <- c("Fixation","Mu")
    temp <- as.data.frame(temp)

    pon <- rbind(pon,temp)
}

pon$Fixation <- as.numeric(pon$Fixation)
pon$qmin <- as.numeric(pon$qmin)
pon$q0 <- as.numeric(pon$q0)
pon$q1 <- as.numeric(pon$q1)
pon$q2 <- as.numeric(pon$q2)
pon$q3 <- as.numeric(pon$q3)
pon$q4 <- as.numeric(pon$q4)
pon$qmax <- as.numeric(pon$qmax)


temp <- pon[which(pon$Fixation > 100),]

image.name <- paste(sep="","Proportion_ON_A",A,"T",T,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- paste(sep="","Proportion ON State Bits From ",100," to ",maxfix+1," Fixations Over ",maxT+1," Trajectories")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",A,", T = ",T,", S = ",s,", n = ",n,", k = ",k,",")
image.subtitles[2] <- paste(sep="","= ",phi,", ",cs,", ",rate,", Initialized ",init,",")

png(image.name,width=2200,height=1200)
print(ggplot(temp,aes(x=Fixation,y=q2))+geom_line(aes(color="black"))+facet_grid(cols=vars(Mu))+geom_ribbon(data=temp,aes(ymin=q1,ymax=q3),color="green",alpha=0.3)+geom_ribbon(data=temp,aes(ymin=q0,ymax=q4),color="purple",alpha=0.2)+geom_ribbon(data=temp,aes(ymin=qmin,ymax=qmax),color="red",alpha=0.1)+geom_hline(yintercept=0.5,color="blue")+theme(text = element_text(size=30),strip.text.y=element_text(size=25),axis.text.y=element_text(size=20),legend.position="bottom")+ylab(bquote("Proportion ON"))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red","y = 0.5"="blue"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema","y = 0.5"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxfix+1)/4),3*floor((maxfix+1)/4),floor((maxfix+1)/4))))
dev.off()


########## Unmasked Proportion for Individuals


library("tidyverse")

codon.frequencies <- c("all","frequencies are 10")
rate <- "rate dep"
initialization <- "all_0s"
init <- "all 0s"
codon.space <- "Cring"
cs <- "codon ring"
s <- 5
phi <- 0.5
k <- 4

MU <- c("0.0001","0.001","0.005","0.01","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2")

Mu <- MU[1]

file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_int.dat")

print(paste(sep="","Loading file ",file))
int0 <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))
print(paste(sep="",file," loaded"))


print(paste(sep="","Ordering ",file))
int0 <- int0[order(int0$Trajectory),]

int0 <- cbind(int0,paste(sep="","Mu = ",rep(Mu,dim(int0)[1])))

colnames(int0)[7] <- "Mu"

int0 <- as_tibble(int0)

n <- nchar(int0$Value[1])
maxfix <- max(int0$Fixation)
maxT <- max(int0$Trajectory)
trnas <- unique(int0$Number[which(int0$Molecule == "tRNA")])
aarss <- unique(int0$Number[which(int0$Molecule == "aaRS")])
ntrna <- length(trnas)
naars <- length(aarss)

print(paste(sep="","Finding decimal form"))
int0 <- int0 %>% mutate(V.dec = strtoi(Value,2))

weight <- function(x){
    return(sum(as.integer(intToBits(x))))
}

int1 <- int0[which(int0$Type == "Mask"),]

print(paste(sep="","Computing weight."))
Weight <- sapply(int1$V.dec,weight)

tempmat <- matrix(Weight,nrow=(maxfix+1)*(ntrna+naars))

print(paste(sep="","Computing quantiles."))
qmin <- apply(tempmat,1,min)
q0 <- apply(tempmat,1,quantile,probs=0.05)
q1 <- apply(tempmat,1,quantile,probs=0.25)
q2 <- apply(tempmat,1,quantile,probs=0.5)
q3 <- apply(tempmat,1,quantile,probs=0.75)
q4 <- apply(tempmat,1,quantile,probs=0.95)
qmax <- apply(tempmat,1,max)

print(paste(sep="","Making Q"))
Q <- cbind(int1[which(int1$Trajectory == 0),c(1:4,7)],qmin,q0,q1,q2,q3,q4,qmax)

for(Mu in MU[-1]){

    file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_int.dat")

    print(paste(sep="","Loading file ",file))
    int0 <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))
    print(paste(sep="",file," loaded"))
    
    int0 <- int0[order(int0$Trajectory),]

    int0 <- cbind(int0,paste(sep="","Mu = ",rep(Mu,dim(int0)[1])))

    colnames(int0)[7] <- "Mu"
    
    int0 <- as_tibble(int0)

    n <- nchar(int0$Value[1])
    maxfix <- max(int0$Fixation)
    maxT <- max(int0$Trajectory)
    trnas <- unique(int0$Number[which(int0$Molecule == "tRNA")])
    aarss <- unique(int0$Number[which(int0$Molecule == "aaRS")])
    ntrna <- length(trnas)
    naars <- length(aarss)

    print(paste(sep="","Finding decimal form"))
    int0 <- int0 %>% mutate(V.dec = strtoi(Value,2))

    int1 <- int0[which(int0$Type == "Mask"),]

    print(paste(sep="","Computing weight."))
    Weight <- sapply(int1$V.dec,weight)

    tempmat <- matrix(Weight,nrow=(maxfix+1)*(ntrna+naars))

    print(paste(sep="","Computing quantiles."))
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    temp <- cbind(int1[which(int1$Trajectory == 0),c(1:4,7)],qmin,q0,q1,q2,q3,q4,qmax)

    print(paste(sep="","Binding"))
    Q <- rbind(Q,temp)
    
}

image.name <- paste(sep="","Proportion_ON_aaRS_A",length(aarss),"T",length(trnas),"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- paste(sep="","Mask Weight for aaRSs Over ",maxT+1," Trajectories")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",naars,", T = ",ntrna,", S = ",s,", n = ",n,", k = ",k,",")
image.subtitles[2] <- paste(sep="","= ",phi,", ",cs,", ",rate,", Initialized ",init,",")

qa <- Q[which(Q$Molecule == "aaRS"),]
qt <- Q[which(Q$Molecule == "tRNA"),]

png(image.name,height=1350,width=2200)
print(ggplot(qa,aes(x=Fixation,y=q2))+facet_grid(rows=vars(Number),cols=vars(Mu))+geom_line(aes(color="black"))+geom_ribbon(data=qa,aes(ymin=q1,ymax=q3),color="green",alpha=0.3)+geom_ribbon(data=qa,aes(ymin=q0,ymax=q4),color="purple",alpha=0.2)+geom_ribbon(data=qa,aes(ymin=qmin,ymax=qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Mask Weight"))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()

image.name <- paste(sep="","Proportion_ON_tRNA_A",length(aarss),"T",length(trnas),"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- paste(sep="","Mask Weight for tRNAs Over ",maxT+1," Trajectories")

png(image.name,height=1350,width=2200)
print(ggplot(qt,aes(x=Fixation,y=q2))+facet_grid(rows=vars(Number),cols=vars(Mu))+geom_line(aes(color="black"))+geom_ribbon(data=qt,aes(ymin=q1,ymax=q3),color="green",alpha=0.3)+geom_ribbon(data=qt,aes(ymin=q0,ymax=q4),color="purple",alpha=0.2)+geom_ribbon(data=qt,aes(ymin=qmin,ymax=qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Mask Weight"))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()

