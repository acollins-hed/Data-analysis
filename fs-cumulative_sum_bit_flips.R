library("tidyverse")

codon.frequencies <- c("all","frequencies are 10")
rate <- "rate dep"
initialization <- "all_0s"
init <- "all 0s"
codon.space <- "Cring"
cs <- "codon ring"
CS <- gsub(" ","_",cs)
s <- 5
phi <- 0.5
k <- 4

MU <- c("0.0001","0.001","0.005","0.01","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2")

file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_10/atINFLATE_A5T5S5n8k4phi0.5Mu",MU[1],"fix1000Cringall0_rate_dep_uniform_aas_int.dat")

int0 <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))

int0 <- int0[order(int0$Trajectory),]

int0 <- as_tibble(int0)

n <- nchar(int0$Value[1])
maxfix <- max(int0$Fixation)
maxT <- max(int0$Trajectory)
trnas <- unique(int0$Number[which(int0$Molecule == "tRNA")])
aarss <- unique(int0$Number[which(int0$Molecule == "aaRS")])
ntrna <- length(trnas)
naars <- length(aarss)

int0 <- int0 %>% mutate(V.dec = strtoi(Value,2))

weight <- function(x){
    return(sum(as.integer(intToBits(x))))
}

State.Flip <- bitwXor(int0$V.dec[which(int0$Fixation > 0 & int0$Type == "State")],int0$V.dec[which(int0$Fixation < maxfix & int0$Type == "State")])
Mask.Flip <- bitwXor(int0$V.dec[which(int0$Fixation > 0 & int0$Type == "Mask")],int0$V.dec[which(int0$Fixation < maxfix & int0$Type == "Mask")])

int2 <- as.data.frame(int0 %>% select(-Value) %>% pivot_wider(names_from=Type,values_from=V.dec))

int2 <- cbind(int2[which(int2$Fixation > 0),],State.Flip,Mask.Flip)

int2 <- cbind(int2,sapply(int2$Mask.Flip,weight))
colnames(int2)[9] <- "Mask.Count"

int2 <- as_tibble(int2)

int2 <- int2 %>% mutate(ON=bitwAnd(Mask,State.Flip))

int2 <- as.data.frame(int2)
int2 <- cbind(int2,sapply(int2$ON,weight))
colnames(int2)[11] <- "S.ON.Flip"

int2 <- as_tibble(int2)
int2 <- int2 %>% mutate(S.OFF.Flip = as.numeric(as.logical(bitwAnd(bitwXor(2^n-1,Mask),State.Flip))))

int2 <- as.data.frame(int2)

cum.ON <- as.vector(apply(matrix(int2$S.ON.Flip,nrow=dim(int2[which(int2$Trajectory == 0),])[1]),2,cumsum))

cum.OFF <- as.vector(apply(matrix(int2$S.OFF.Flip,nrow=dim(int2[which(int2$Trajectory == 0),])[1]),2,cumsum))

cum.mask <- as.vector(apply(matrix(int2$Mask.Count,nrow=dim(int2[which(int2$Trajectory == 0),])[1]),2,cumsum))

int2 <- cbind(int2,rep(paste(sep="","Mu = ",MU[1]),dim(int2)[1]))
colnames(int2)[dim(int2)[2]] <- "Mu"

int3 <- cbind(int2[,1:6],cum.ON,cum.OFF,cum.mask)

int <- int3 %>% mutate(state.diff=cum.OFF-cum.ON,ON.diff = cum.ON-cum.mask, OFF.diff = cum.OFF-cum.mask)

tempmat <- matrix(int$ON.diff,nrow=dim(int[which(int$Trajectory == 0),])[1])
ON.qmin <- apply(tempmat,1,min)
ON.q0 <- apply(tempmat,1,quantile,probs=0.05)
ON.q1 <- apply(tempmat,1,quantile,probs=0.25)
ON.q2 <- apply(tempmat,1,quantile,probs=0.5)
ON.q3 <- apply(tempmat,1,quantile,probs=0.75)
ON.q4 <- apply(tempmat,1,quantile,probs=0.95)
ON.qmax <- apply(tempmat,1,max)

tempmat <- matrix(int$OFF.diff,nrow=dim(int[which(int$Trajectory == 0),])[1])
OFF.qmin <- apply(tempmat,1,min)
OFF.q0 <- apply(tempmat,1,quantile,probs=0.05)
OFF.q1 <- apply(tempmat,1,quantile,probs=0.25)
OFF.q2 <- apply(tempmat,1,quantile,probs=0.5)
OFF.q3 <- apply(tempmat,1,quantile,probs=0.75)
OFF.q4 <- apply(tempmat,1,quantile,probs=0.95)
OFF.qmax <- apply(tempmat,1,max)

tempmat <- matrix(int$state.diff,nrow=dim(int[which(int$Trajectory == 0),])[1])
state.qmin <- apply(tempmat,1,min)
state.q0 <- apply(tempmat,1,quantile,probs=0.05)
state.q1 <- apply(tempmat,1,quantile,probs=0.25)
state.q2 <- apply(tempmat,1,quantile,probs=0.5)
state.q3 <- apply(tempmat,1,quantile,probs=0.75)
state.q4 <- apply(tempmat,1,quantile,probs=0.95)
state.qmax <- apply(tempmat,1,max)

int1 <- cbind(int[which(int$Trajectory == 0),1:6],ON.qmin,ON.q0,ON.q1,ON.q2,ON.q3,ON.q4,ON.qmax,OFF.qmin,OFF.q0,OFF.q1,OFF.q2,OFF.q3,OFF.q4,OFF.qmax,state.qmin,state.q0,state.q1,state.q2,state.q3,state.q4,state.qmax,rep(paste(sep="","Mu = ",MU[1]),length(ON.qmin)))

colnames(int1)[28] <- "Mu"

for(M in MU[-1]){

    file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",M,"fix1000C2all0_rate_dep_uniform_aas_int.dat")

    int0 <- read.table(file,header=TRUE,row.names=NULL,colClasses=c("numeric","numeric","character","character","character","character"))

    int0 <- int0[order(int0$Trajectory),]

    int0 <- as_tibble(int0)
    
    n <- nchar(int0$Value[1])
    maxfix <- max(int0$Fixation)
    maxT <- max(int0$Trajectory)
    trnas <- unique(int0$Number[which(int0$Molecule == "tRNA")])
    aarss <- unique(int0$Number[which(int0$Molecule == "aaRS")])
    ntrna <- length(trnas)
    naars <- length(aarss)
    
    int0 <- int0 %>% mutate(V.dec = strtoi(Value,2))
    
    State.Flip <- bitwXor(int0$V.dec[which(int0$Fixation > 0 & int0$Type == "State")],int0$V.dec[which(int0$Fixation < maxfix & int0$Type == "State")])
Mask.Flip <- bitwXor(int0$V.dec[which(int0$Fixation > 0 & int0$Type == "Mask")],int0$V.dec[which(int0$Fixation < maxfix & int0$Type == "Mask")])
2
    temp2 <- as.data.frame(int0 %>% select(-Value) %>% pivot_wider(names_from=Type,values_from=V.dec))
    
    temp2 <- cbind(temp2[which(temp2$Fixation > 0),],State.Flip,Mask.Flip)
    
    temp2 <- cbind(temp2,sapply(temp2$Mask.Flip,weight))
    colnames(temp2)[9] <- "Mask.Count"
    
    temp2 <- as_tibble(temp2)
    
    temp2 <- temp2 %>% mutate(ON=bitwAnd(Mask,State.Flip))
    
    temp2 <- as.data.frame(temp2)
    temp2 <- cbind(temp2,sapply(temp2$ON,weight))
    colnames(temp2)[11] <- "S.ON.Flip"
    
    temp2 <- as_tibble(temp2)
    temp2 <- temp2 %>% mutate(S.OFF.Flip = as.numeric(as.logical(bitwAnd(bitwXor(2^n-1,Mask),State.Flip))))
    
    temp2 <- as.data.frame(temp2)
    
    cum.ON <- as.vector(apply(matrix(temp2$S.ON.Flip,nrow=dim(temp2[which(temp2$Trajectory == 0),])[1]),2,cumsum))
    
    cum.OFF <- as.vector(apply(matrix(temp2$S.OFF.Flip,nrow=dim(temp2[which(temp2$Trajectory == 0),])[1]),2,cumsum))
    
    cum.mask <- as.vector(apply(matrix(temp2$Mask.Count,nrow=dim(temp2[which(temp2$Trajectory == 0),])[1]),2,cumsum))
    
    int3 <- cbind(temp2[,1:6],cum.ON,cum.OFF,cum.mask)

    temp2 <- cbind(temp2,rep(paste(sep="","Mu = ",M),dim(temp2)[1]))
    colnames(temp2)[dim(temp2)[2]] <- "Mu"

    int2 <- rbind(int2,temp2)
    
    int <- int3 %>% mutate(state.diff=cum.OFF-cum.ON,ON.diff = cum.ON-cum.mask, OFF.diff = cum.OFF-cum.mask)
    
    tempmat <- matrix(int$ON.diff,nrow=dim(int[which(int$Trajectory == 0),])[1])
    ON.qmin <- apply(tempmat,1,min)
    ON.q0 <- apply(tempmat,1,quantile,probs=0.05)
    ON.q1 <- apply(tempmat,1,quantile,probs=0.25)
    ON.q2 <- apply(tempmat,1,quantile,probs=0.5)
    ON.q3 <- apply(tempmat,1,quantile,probs=0.75)
    ON.q4 <- apply(tempmat,1,quantile,probs=0.95)
    ON.qmax <- apply(tempmat,1,max)
    
    tempmat <- matrix(int$OFF.diff,nrow=dim(int[which(int$Trajectory == 0),])[1])
    OFF.qmin <- apply(tempmat,1,min)
    OFF.q0 <- apply(tempmat,1,quantile,probs=0.05)
    OFF.q1 <- apply(tempmat,1,quantile,probs=0.25)
    OFF.q2 <- apply(tempmat,1,quantile,probs=0.5)
    OFF.q3 <- apply(tempmat,1,quantile,probs=0.75)
    OFF.q4 <- apply(tempmat,1,quantile,probs=0.95)
    OFF.qmax <- apply(tempmat,1,max)
    
    tempmat <- matrix(int$state.diff,nrow=dim(int[which(int$Trajectory == 0),])[1])
    state.qmin <- apply(tempmat,1,min)
    state.q0 <- apply(tempmat,1,quantile,probs=0.05)
    state.q1 <- apply(tempmat,1,quantile,probs=0.25)
    state.q2 <- apply(tempmat,1,quantile,probs=0.5)
    state.q3 <- apply(tempmat,1,quantile,probs=0.75)
    state.q4 <- apply(tempmat,1,quantile,probs=0.95)
    state.qmax <- apply(tempmat,1,max)

    temp <- cbind(int[which(int$Trajectory == 0),1:6],ON.qmin,ON.q0,ON.q1,ON.q2,ON.q3,ON.q4,ON.qmax,OFF.qmin,OFF.q0,OFF.q1,OFF.q2,OFF.q3,OFF.q4,OFF.qmax,state.qmin,state.q0,state.q1,state.q2,state.q3,state.q4,state.qmax,rep(paste(sep="","Mu = ",M),length(ON.qmin)))

    colnames(temp)[28] <- "Mu"

    int1 <- rbind(int1,temp)
}

image.name <- paste(sep="","Masked_State-Mask_A",length(aarss),"T",length(trnas),"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- paste(sep="","Difference of Cumulative Sum of Masked State Bit Fixations and Mask Bit Fixations Over ",maxT+1," Trajectories")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",naars,", T = ",ntrna,", S = ",s,", n = ",n,", k = ",k,",")
image.subtitles[2] <- paste(sep="","= ",phi,", ",cs,", ",rate,", Initialized ",init,",")

png(image.name,height=1200,width=2200)
print(ggplot(int1,aes(x=Fixation,y=OFF.q2))+facet_wrap(~Mu)+geom_line(aes(color="black"))+geom_ribbon(data=int1,aes(ymin=OFF.q1,ymax=OFF.q3),color="green",alpha=0.3)+geom_ribbon(data=int1,aes(ymin=OFF.q0,ymax=OFF.q4),color="purple",alpha=0.2)+geom_ribbon(data=int1,aes(ymin=OFF.qmin,ymax=OFF.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Masked"~Delta))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2]))))
dev.off()

image.name <- paste(sep="","Masked_State-Mask_A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,"_alt.png")

png(image.name,height=1200,width=2200)
print(ggplot(int1,aes(x=Fixation,y=OFF.q2))+facet_wrap(~Mu,nrow=1)+geom_line(aes(color="black"))+geom_ribbon(data=int1,aes(ymin=OFF.q1,ymax=OFF.q3),color="green",alpha=0.3)+geom_ribbon(data=int1,aes(ymin=OFF.q0,ymax=OFF.q4),color="purple",alpha=0.2)+geom_ribbon(data=int1,aes(ymin=OFF.qmin,ymax=OFF.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Masked"~Delta))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2]))))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4)))
dev.off()

image.name <- paste(sep="","Unmasked_State-Mask_A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- paste(sep="","Difference of Cumulative Sum of Unmasked State Bit Fixations and Mask Bit Fixations Over ",maxT+1," Trajectories")

png(image.name,height=1200,width=2200)
print(ggplot(int1,aes(x=Fixation,y=ON.q2))+facet_wrap(~Mu)+geom_line(aes(color="black"))+geom_ribbon(data=int1,aes(ymin=ON.q1,ymax=ON.q3),color="green",alpha=0.3)+geom_ribbon(data=int1,aes(ymin=ON.q0,ymax=ON.q4),color="purple",alpha=0.2)+geom_ribbon(data=int1,aes(ymin=ON.qmin,ymax=ON.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Unmasked"~Delta))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2]))))
dev.off()

image.name <- paste(sep="","Unmasked_State-Mask_A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,"_alt.png")
png(image.name,height=1200,width=2200)
print(ggplot(int1,aes(x=Fixation,y=ON.q2))+facet_wrap(~Mu,nrow=1)+geom_line(aes(color="black"))+geom_ribbon(data=int1,aes(ymin=ON.q1,ymax=ON.q3),color="green",alpha=0.3)+geom_ribbon(data=int1,aes(ymin=ON.q0,ymax=ON.q4),color="purple",alpha=0.2)+geom_ribbon(data=int1,aes(ymin=ON.qmin,ymax=ON.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Unmasked"~Delta))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()

image.name <- paste(sep="","Masked-Unmasked_State_A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- paste(sep="","Difference of Cumulative Sum of Masked State Bit Fixations and Unmasked State Bit Fixations and Mask Bit Fixations Over ",maxT+1," Trajectories")


png(image.name,height=1200,width=2200)
print(ggplot(int1,aes(x=Fixation,y=state.q2))+facet_wrap(~Mu)+geom_line(aes(color="black"))+geom_ribbon(data=int1,aes(ymin=state.q1,ymax=state.q3),color="green",alpha=0.3)+geom_ribbon(data=int1,aes(ymin=state.q0,ymax=state.q4),color="purple",alpha=0.2)+geom_ribbon(data=int1,aes(ymin=state.qmin,ymax=state.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote({Delta}[Masked]-{Delta}[Unmaksed]))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2]))))
dev.off()

image.name <- paste(sep="","Masked-Unmasked_State_A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,"_alt.png")

png(image.name,height=1200,width=2200)
print(ggplot(int1,aes(x=Fixation,y=state.q2))+facet_wrap(~Mu,nrow=1)+geom_line(aes(color="black"))+geom_ribbon(data=int1,aes(ymin=state.q1,ymax=state.q3),color="green",alpha=0.3)+geom_ribbon(data=int1,aes(ymin=state.q0,ymax=state.q4),color="purple",alpha=0.2)+geom_ribbon(data=int1,aes(ymin=state.qmin,ymax=state.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote({Delta}[Masked]-{Delta}[Unmaksed]))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()


int4 <- cbind(int2[which(int2$Molecule == "tRNA" & int2$Number == trnas[1]),1:6],as.vector(apply(matrix(int2$S.ON.Flip[which(int2$Molecule == "tRNA" & int2$Number == trnas[1])],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == trnas[1] & int2$Molecule == "tRNA" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),as.vector(apply(matrix(int2$S.OFF.Flip[which(int2$Molecule == "tRNA" & int2$Number == trnas[1])],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == trnas[1] & int2$Molecule == "tRNA" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),as.vector(apply(matrix(int2$Mask.Count[which(int2$Molecule == "tRNA" & int2$Number == trnas[1])],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == trnas[1] & int2$Molecule == "tRNA" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),int2$Mu[which(int2$Molecule == "tRNA" & int2$Number == trnas[1])])
colnames(int4)[7:10] <- c("cum.ON","cum.OFF","cum.mask","Mu")

for(t in trnas[-1]){

    temp <- cbind(int2[which(int2$Molecule == "tRNA" & int2$Number == t),1:6],as.vector(apply(matrix(int2$S.ON.Flip[which(int2$Molecule == "tRNA" & int2$Number == t)],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == t & int2$Molecule == "tRNA" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),as.vector(apply(matrix(int2$S.OFF.Flip[which(int2$Molecule == "tRNA" & int2$Number == t)],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == t & int2$Molecule == "tRNA" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),as.vector(apply(matrix(int2$Mask.Count[which(int2$Molecule == "tRNA" & int2$Number == t)],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == t & int2$Molecule == "tRNA" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),int2$Mu[which(int2$Molecule == "tRNA" & int2$Number == t)])

    colnames(temp)[7:10] <- c("cum.ON","cum.OFF","cum.mask","Mu")
    
    int4 <- rbind(int4,temp)
}

temp <- cbind(int2[which(int2$Molecule == "aaRS" & int2$Number == aarss[1]),1:6],as.vector(apply(matrix(int2$S.ON.Flip[which(int2$Molecule == "aaRS" & int2$Number == aarss[1])],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == aarss[1] & int2$Molecule == "aaRS" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),as.vector(apply(matrix(int2$S.OFF.Flip[which(int2$Molecule == "aaRS" & int2$Number == aarss[1])],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == aarss[1] & int2$Molecule == "aaRS" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),as.vector(apply(matrix(int2$Mask.Count[which(int2$Molecule == "aaRS" & int2$Number == aarss[1])],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == aarss[1] & int2$Molecule == "aaRS" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),int2$Mu[which(int2$Molecule == "aaRS" & int2$Number == aarss[1])])

colnames(temp)[7:10] <- c("cum.ON","cum.OFF","cum.mask","Mu")


int4 <- rbind(int4,temp)

for(a in aarss[-1]){

    temp <- cbind(int2[which(int2$Molecule == "aaRS" & int2$Number == a),1:6],as.vector(apply(matrix(int2$S.ON.Flip[which(int2$Molecule == "aaRS" & int2$Number == a)],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == a & int2$Molecule == "aaRS" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),as.vector(apply(matrix(int2$S.OFF.Flip[which(int2$Molecule == "aaRS" & int2$Number == a)],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == a & int2$Molecule == "aaRS" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),as.vector(apply(matrix(int2$Mask.Count[which(int2$Molecule == "aaRS" & int2$Number == a)],nrow=dim(int2[which(int2$Trajectory == 0 & int2$Number == a & int2$Molecule == "aaRS" & int2$Mu == paste(sep="","Mu = ",MU[1])),])[1]),2,cumsum)),int2$Mu[which(int2$Molecule == "aaRS" & int2$Number == a)])

    colnames(temp)[7:10] <- c("cum.ON","cum.OFF","cum.mask","Mu")
    
    int4 <- rbind(int4,temp)
}

int4 <- as_tibble(int4)
int4 <- int4 %>% mutate(state.diff = cum.OFF - cum.ON,OFF.diff = cum.OFF-cum.mask,ON.diff=cum.ON-cum.mask)

tempmat <- matrix(int4$ON.diff,nrow=dim(int4[which(int4$Trajectory == 0),])[1])

ON.qmin <- apply(tempmat,1,min)
ON.q0 <- apply(tempmat,1,quantile,probs=0.05)
ON.q1 <- apply(tempmat,1,quantile,probs=0.25)
ON.q2 <- apply(tempmat,1,quantile,probs=0.5)
ON.q3 <- apply(tempmat,1,quantile,probs=0.75)
ON.q4 <- apply(tempmat,1,quantile,probs=0.95)
ON.qmax <- apply(tempmat,1,max)

tempmat <- matrix(int4$OFF.diff,nrow=dim(int4[which(int4$Trajectory == 0),])[1])

OFF.qmin <- apply(tempmat,1,min)
OFF.q0 <- apply(tempmat,1,quantile,probs=0.05)
OFF.q1 <- apply(tempmat,1,quantile,probs=0.25)
OFF.q2 <- apply(tempmat,1,quantile,probs=0.5)
OFF.q3 <- apply(tempmat,1,quantile,probs=0.75)
OFF.q4 <- apply(tempmat,1,quantile,probs=0.95)
OFF.qmax <- apply(tempmat,1,max)

tempmat <- matrix(int4$state.diff,nrow=dim(int4[which(int4$Trajectory == 0),])[1])

state.qmin <- apply(tempmat,1,min)
state.q0 <- apply(tempmat,1,quantile,probs=0.05)
state.q1 <- apply(tempmat,1,quantile,probs=0.25)
state.q2 <- apply(tempmat,1,quantile,probs=0.5)
state.q3 <- apply(tempmat,1,quantile,probs=0.75)
state.q4 <- apply(tempmat,1,quantile,probs=0.95)
state.qmax <- apply(tempmat,1,max)


int5 <- cbind(int4[which(int4$Trajectory == 0),1:6],ON.qmin,ON.q0,ON.q1,ON.q2,ON.q3,ON.q4,ON.qmax,OFF.qmin,OFF.q0,OFF.q1,OFF.q2,OFF.q3,OFF.q4,OFF.qmax,state.qmin,state.q0,state.q1,state.q2,state.q3,state.q4,state.qmax,int4$Mu[which(int4$Trajectory == 0)])

colnames(int5)[dim(int5)[2]] <- "Mu"

image.name <- paste(sep="","Masked_State-Mask_aaRSs-A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,"_alt.png")
image.title <- paste(sep="","Difference of Cumulative Sum of aaRS Masked State Bit Fixations and Mask Bit Fixations Over ",maxT+1," Trajectories")

png(image.name,height=1200,width=2200)
print(ggplot(int5[which(int5$Molecule == "aaRS"),],aes(x=Fixation,y=OFF.q2))+geom_line(aes(color="black"))+facet_grid(cols=vars(Mu),rows=vars(Number))+geom_ribbon(data=int5[which(int5$Molecule == "aaRS"),],aes(ymin=OFF.q1,ymax=OFF.q3),color="green",alpha=0.3)+geom_ribbon(data=int5[which(int5$Molecule == "aaRS"),],aes(ymin=OFF.q0,ymax=OFF.q4),color="purple",alpha=0.2)+geom_ribbon(data=int5[which(int5$Molecule == "aaRS"),],aes(ymin=OFF.qmin,ymax=OFF.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Masked"~Delta^"aaRS"))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()

image.name <- paste(sep="","Masked_State-Mask_tRNAs-A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,"_alt.png")
image.title <- paste(sep="","Difference of Cumulative Sum of tRNA Masked State Bit Fixations and Mask Bit Fixations Over ",maxT+1," Trajectories")

png(image.name,height=1200,width=2200)
print(ggplot(int5[which(int5$Molecule == "tRNA"),],aes(x=Fixation,y=OFF.q2))+geom_line(aes(color="black"))+facet_grid(cols=vars(Mu),rows=vars(Number))+geom_ribbon(data=int5[which(int5$Molecule == "tRNA"),],aes(ymin=OFF.q1,ymax=OFF.q3),color="green",alpha=0.3)+geom_ribbon(data=int5[which(int5$Molecule == "tRNA"),],aes(ymin=OFF.q0,ymax=OFF.q4),color="purple",alpha=0.2)+geom_ribbon(data=int5[which(int5$Molecule == "tRNA"),],aes(ymin=OFF.qmin,ymax=OFF.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Masked"~Delta^"tRNA"))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()



image.name <- paste(sep="","Unmasked_State-Mask_aaRSs-A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,"_alt.png")
image.title <- paste(sep="","Difference of Cumulative Sum of aaRS Unmasked State Bit Fixations and Mask Bit Fixations Over ",maxT+1," Trajectories")

png(image.name,height=1200,width=2200)
print(ggplot(int5[which(int5$Molecule == "aaRS"),],aes(x=Fixation,y=ON.q2))+geom_line(aes(color="black"))+facet_grid(cols=vars(Mu),rows=vars(Number))+geom_ribbon(data=int5[which(int5$Molecule == "aaRS"),],aes(ymin=ON.q1,ymax=ON.q3),color="green",alpha=0.3)+geom_ribbon(data=int5[which(int5$Molecule == "aaRS"),],aes(ymin=ON.q0,ymax=ON.q4),color="purple",alpha=0.2)+geom_ribbon(data=int5[which(int5$Molecule == "aaRS"),],aes(ymin=ON.qmin,ymax=ON.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Unmasked"~Delta^"aaRS"))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()

image.name <- paste(sep="","Unmasked_State-Mask_tRNAs-A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,"_alt.png")
image.title <- paste(sep="","Difference of Cumulative Sum of tRNA Unmasked State Bit Fixations and Mask Bit Fixations Over ",maxT+1," Trajectories")

png(image.name,height=1200,width=2200)
print(ggplot(int5[which(int5$Molecule == "tRNA"),],aes(x=Fixation,y=ON.q2))+geom_line(aes(color="black"))+facet_grid(cols=vars(Mu),rows=vars(Number))+geom_ribbon(data=int5[which(int5$Molecule == "tRNA"),],aes(ymin=ON.q1,ymax=ON.q3),color="green",alpha=0.3)+geom_ribbon(data=int5[which(int5$Molecule == "tRNA"),],aes(ymin=ON.q0,ymax=ON.q4),color="purple",alpha=0.2)+geom_ribbon(data=int5[which(int5$Molecule == "tRNA"),],aes(ymin=ON.qmin,ymax=ON.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote("Unmasked"~Delta^"tRNA"))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()


image.name <- paste(sep="","Masked-Unmasked_aaRSs-A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,"_alt.png")
image.title <- paste(sep="","Difference of Cumulative Sum of aaRS Masked State Bit Fixations and Unmasked State Bit Fixations Over ",maxT+1," Trajectories")

png(image.name,height=1200,width=2200)
print(ggplot(int5[which(int5$Molecule == "aaRS"),],aes(x=Fixation,y=state.q2))+geom_line(aes(color="black"))+facet_grid(cols=vars(Mu),rows=vars(Number))+geom_ribbon(data=int5[which(int5$Molecule == "aaRS"),],aes(ymin=state.q1,ymax=state.q3),color="green",alpha=0.3)+geom_ribbon(data=int5[which(int5$Molecule == "aaRS"),],aes(ymin=state.q0,ymax=state.q4),color="purple",alpha=0.2)+geom_ribbon(data=int5[which(int5$Molecule == "aaRS"),],aes(ymin=state.qmin,ymax=state.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote({Delta^"aaRS"}[Masked]-{Delta^"aaRS"}[Unmasked]))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxT+1)/4),3*floor((maxT+1)/4),floor((maxT+1)/4))))
dev.off()

image.name <- paste(sep="","Masked-Unmasked_tRNAs-A",naars,"T",ntrna,"S",s,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,"_alt.png")
image.title <- paste(sep="","Difference of Cumulative Sum of tRNA Masked State Bit Fixations and Unmasked State Bit Fixations Over ",maxT+1," Trajectories")

png(image.name,height=1200,width=2200)
print(ggplot(int5[which(int5$Molecule == "tRNA"),],aes(x=Fixation,y=state.q2))+geom_line(aes(color="black"))+facet_grid(cols=vars(Mu),rows=vars(Number))+geom_ribbon(data=int5[which(int5$Molecule == "tRNA"),],aes(ymin=state.q1,ymax=state.q3),color="green",alpha=0.3)+geom_ribbon(data=int5[which(int5$Molecule == "tRNA"),],aes(ymin=state.q0,ymax=state.q4),color="purple",alpha=0.2)+geom_ribbon(data=int5[which(int5$Molecule == "tRNA"),],aes(ymin=state.qmin,ymax=state.qmax),color="red",alpha=0.1)+theme(text = element_text(size=30),legend.position="bottom")+ylab(bquote({Delta^"tRNA"}[Masked]-{Delta^"tRNA"}[Unmasked]))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_x_continuous(breaks=seq(floor((maxfix+1)/4),3*floor((maxfix+1)/4),floor((maxfix+1)/4))))
dev.off()

