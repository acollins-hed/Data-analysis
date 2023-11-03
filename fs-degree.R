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
n <- 8

MU <- c("0.0001","0.001","0.005","0.01","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2")

Mu <- MU[1]

file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_code.dat")

print(paste(sep="","Loading file ",file))
code <- read.table(file,header=TRUE,row.names=NULL)
print(paste(sep="",file," loaded"))

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
temp <- code %>% select(-Prob_Interaction) %>% pivot_wider(names_from=tRNA,values_from=Match)

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
Qa <- cbind(rep(0:maxfix,each=naars),rep(aarss,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","Mu = ",Mu),length(q2)))

Qa <- as.data.frame(Qa)

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
Qt <- cbind(rep(0:maxfix,each=ntrna),rep(trnas,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","Mu = ",Mu),length(q2)))

Qt <- as.data.frame(Qt)

for(Mu in MU[-1]){

    file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_code.dat")

    print(paste(sep="","Loading file ",file))
    code <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," loaded"))
    
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
    temp <- code %>% select(-Prob_Interaction) %>% pivot_wider(names_from=tRNA,values_from=Match)
    
    print(paste(sep="","Finding each aaRS degree."))
    adegree <- apply(temp[,-1:-3],1,sum)
    
    tempmat <- matrix(adegree,ncol=maxT+1)
    
    print(paste(sep="","Finding quantiles."))
    qmin <- apply(tempmat,1,min)
    q0 <- apply(tempmat,1,quantile,probs=0.05)
    q1 <- apply(tempmat,1,quantile,probs=0.25)
    q2 <- apply(tempmat,1,quantile,probs=0.5)
    q3 <- apply(tempmat,1,quantile,probs=0.75)
    q4 <- apply(tempmat,1,quantile,probs=0.95)
    qmax <- apply(tempmat,1,max)
    
    print(paste(sep="","Making Qa."))
    temp <- cbind(rep(0:maxfix,each=naars),rep(aarss,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","Mu = ",Mu),length(q2)))
    
    temp <- as.data.frame(temp)
    
    Qa <- rbind(Qa,temp)

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
    temp <- cbind(rep(0:maxfix,each=ntrna),rep(trnas,maxfix+1),qmin,q0,q1,q2,q3,q4,qmax,rep(paste(sep="","Mu = ",Mu),length(q2)))
    
    temp <- as.data.frame(temp)

    Qt <- rbind(Qt,temp)
}

colnames(Qa)[c(1,2,10)] <- c("Fixation","aaRS","Mu")

Qa$Fixation <- as.numeric(Qa$Fixation)
Qa$qmin <- as.numeric(Qa$qmin)
Qa$q0 <- as.numeric(Qa$q0)
Qa$q1 <- as.numeric(Qa$q1)
Qa$q2 <- as.numeric(Qa$q2)
Qa$q3 <- as.numeric(Qa$q3)
Qa$q4 <- as.numeric(Qa$q4)
Qa$qmax <- as.numeric(Qa$qmax)


colnames(Qt)[c(1,2,10)] <- c("Fixation","tRNA","Mu")

Qt$Fixation <- as.numeric(Qt$Fixation)
Qt$qmin <- as.numeric(Qt$qmin)
Qt$q0 <- as.numeric(Qt$q0)
Qt$q1 <- as.numeric(Qt$q1)
Qt$q2 <- as.numeric(Qt$q2)
Qt$q3 <- as.numeric(Qt$q3)
Qt$q4 <- as.numeric(Qt$q4)
Qt$qmax <- as.numeric(Qt$qmax)

image.name <- paste(sep="","Degree_aaRSs_A",length(aarss),"T",length(trnas),"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- paste(sep="","Degree Quantiles for each aaRS Over ",maxT+1," Trajectories, ",maxfix+1," Fixations Each")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",naars,", T = ",ntrna,", S = ",s,", n = ",n,", k = ",k,",")
image.subtitles[2] <- paste(sep="","= ",phi,", ",cs,", ",rate,", Initialized ",init,",")

png(image.name,height=1350,width=2200)
print(ggplot(Qa,aes(x=Fixation,y=q2))+facet_grid(rows=vars(aaRS),cols=vars(Mu))+geom_line(aes(color="black"))+geom_ribbon(data=Qa,aes(ymin=q1,ymax=q3),color="green",alpha=0.3)+geom_ribbon(data=Qa,aes(ymin=q0,ymax=q4),color="purple",alpha=0.2)+geom_ribbon(data=Qa,aes(ymin=qmin,ymax=qmax),color="red",alpha=0.1)+theme(text = element_text(size=35),axis.text=element_text(size=20),legend.position="bottom")+ylab(bquote("Degree"))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red","y = 0.5"="blue"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema","y = 4"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxfix+1)/4),3*floor((maxfix+1)/4),floor((maxfix+1)/4)))+geom_hline(yintercept=4,color="blue"))
dev.off()


image.name <- paste(sep="","Degree_tRNAs_A",length(aarss),"T",length(trnas),"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- paste(sep="","Degree Quantiles for each tRNA Over ",maxT+1," Trajectories, ",maxfix+1," Fixations Each")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",naars,", T = ",ntrna,", S = ",s,", n = ",n,", k = ",k,",")
image.subtitles[2] <- paste(sep="","= ",phi,", ",cs,", ",rate,", Initialized ",init,",")

png(image.name,height=1350,width=2200)
print(ggplot(Qt,aes(x=Fixation,y=q2))+facet_grid(rows=vars(tRNA),cols=vars(Mu))+geom_line(aes(color="black"))+geom_ribbon(data=Qt,aes(ymin=q1,ymax=q3),color="green",alpha=0.3)+geom_ribbon(data=Qt,aes(ymin=q0,ymax=q4),color="purple",alpha=0.2)+geom_ribbon(data=Qt,aes(ymin=qmin,ymax=qmax),color="red",alpha=0.1)+theme(text = element_text(size=35),axis.text=element_text(size=20),legend.position="bottom")+ylab(bquote("Degree"))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red","y = 0.5"="blue"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema","y = 4"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxfix+1)/4),3*floor((maxfix+1)/4),floor((maxfix+1)/4)))+geom_hline(yintercept=4,color="blue"))
dev.off()

