codon.frequencies <- c("all","frequencies are 10")
rate <- "rate dep"
Rate <- gsub(" ","_",rate)
initialization <- "all_0s"
init <- gsub("_"," ",initialization)
init.alt <- gsub("s","",gsub(" ","",init)) 
Codon.space <- "Cring"
codon.space <- "cring"
cs <- "codon ring"
phi <- 0.5
k <- 4
n <- 8
naars <- 5
ntrna <- 5
sts <- 5
directory <- "/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/"


MU <- c("0.0001","0.001","0.005","0.01","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2")
#MU <- c("0.001","0.01","0.1","0.2")

Mu <- MU[1]

file <- paste(sep="",directory,"atINFLATE_A",naars,"T",ntrna,"S",sts,"n",n,"k",k,"phi",phi,"Mu",Mu,"fix1000",Codon.space,init.alt,"_",Rate,"_unif_aas_prob.dat")

print(paste(sep="","Loading ",file))
prob <- read.table(file,header=TRUE,row.names=NULL)
print(paste(sep="",file," Loaded"))

AAs <- unique(prob$Amino_Acid)
aas <- length(AAs)
Sts <- unique(prob$Site_Type)
sts <- length(Sts)
p.of.s <- 1/sts

library(ggplot2)

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

print(paste(sep="","Beginning calculations for amino acid ",AAs[1]))
mat <- matrix(p.of.a$Probability[which(p.of.a$Amino_Acid == AAs[1])],nrow=maxfix+1)
qmin <- apply(mat,1,min)
q0 <- apply(mat,1,quantile,probs=0.05)
q1 <- apply(mat,1,quantile,probs=0.25)
q2 <- apply(mat,1,quantile,probs=0.5)
q3 <- apply(mat,1,quantile,probs=0.75)
q4 <- apply(mat,1,quantile,probs=0.95)
qmax <- apply(mat,1,max)

q <- cbind(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(AAs[1],maxfix+1),rep(paste(sep="","Mu = ",Mu),length(qmin)))
q <- as.data.frame(q)
colnames(q)[1] <- "Fixation"
colnames(q)[9] <- "Amino_Acid"
colnames(q)[10] <- "Mu"

print(paste(sep="","Beginning calculations for remaining amino acids"))
for(a in 2:aas){
    mat <- matrix(p.of.a$Probability[which(p.of.a$Amino_Acid == AAs[a])],nrow=maxfix+1)
    qmin <- apply(mat,1,min)
    q0 <- apply(mat,1,quantile,probs=0.05)
    q1 <- apply(mat,1,quantile,probs=0.25)
    q2 <- apply(mat,1,quantile,probs=0.5)
    q3 <- apply(mat,1,quantile,probs=0.75)
    q4 <- apply(mat,1,quantile,probs=0.95)
    qmax <- apply(mat,1,max)

    temp <- cbind(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(AAs[a],maxfix+1),rep(paste(sep="","Mu = ",Mu),length(qmin)))
    temp <- as.data.frame(temp)
    colnames(temp)[1] <- "Fixation"
    colnames(temp)[9] <- "Amino_Acid"
    colnames(temp)[10] <- "Mu"
    q <- rbind(q,temp)
}


for(Mu in MU[-1]){
    
    file <- paste(sep="",directory,"atINFLATE_A",naars,"T",ntrna,"S",sts,"n",n,"k",k,"phi",phi,"Mu",Mu,"fix1000",Codon.space,init.alt,"_",Rate,"_unif_aas_prob.dat")

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

    print(paste(sep="","Beginning calculations for amino acid ",AAs[1]))
    mat <- matrix(p.of.a$Probability[which(p.of.a$Amino_Acid == AAs[1])],nrow=maxfix+1)
    qmin <- apply(mat,1,min)
    q0 <- apply(mat,1,quantile,probs=0.05)
    q1 <- apply(mat,1,quantile,probs=0.25)
    q2 <- apply(mat,1,quantile,probs=0.5)
    q3 <- apply(mat,1,quantile,probs=0.75)
    q4 <- apply(mat,1,quantile,probs=0.95)
    qmax <- apply(mat,1,max)
    
    temp <- cbind(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(AAs[1],maxfix+1),rep(paste(sep="","Mu = ",Mu),length(qmin)))
    temp <- as.data.frame(temp)
    colnames(temp)[1] <- "Fixation"
    colnames(temp)[9] <- "Amino_Acid"
    colnames(temp)[10] <- "Mu"

    print(paste(sep="","Beginning calculations for remaining amino acids"))
    for(a in 2:aas){
        mat <- matrix(p.of.a$Probability[which(p.of.a$Amino_Acid == AAs[a])],nrow=maxfix+1)
        qmin <- apply(mat,1,min)
        q0 <- apply(mat,1,quantile,probs=0.05)
        q1 <- apply(mat,1,quantile,probs=0.25)
        q2 <- apply(mat,1,quantile,probs=0.5)
        q3 <- apply(mat,1,quantile,probs=0.75)
        q4 <- apply(mat,1,quantile,probs=0.95)
        qmax <- apply(mat,1,max)
        
        temp1 <- cbind(0:maxfix,qmin,q0,q1,q2,q3,q4,qmax,rep(AAs[a],maxfix+1),rep(paste(sep="","Mu = ",Mu),length(qmin)))
        temp1 <- as.data.frame(temp1)
        colnames(temp1)[1] <- "Fixation"
        colnames(temp1)[9] <- "Amino_Acid"
        colnames(temp1)[10] <- "Mu"
        temp <- rbind(temp,temp1)
    }
    print(paste(sep="","Finishing up ",file))
    q <- rbind(q,temp)
}

q$Fixation <- as.numeric(q$Fixation)
q$qmin <- as.numeric(q$qmin)
q$q0 <- as.numeric(q$q0)
q$q1 <- as.numeric(q$q1)
q$q2 <- as.numeric(q$q2)
q$q3 <- as.numeric(q$q3)
q$q4 <- as.numeric(q$q4)
q$qmax <- as.numeric(q$qmax)

print("Writing q")
ofile <- paste(sep="","A",naars,"T",ntrna,"S",sts,"_phi",phi,"_",init.alt,"_",Rate,"_",codon.space,"_P_Of_AA_Quantiles.dat")
write.table(q,file=ofile,row.names=FALSE)

image.name <- paste(sep="","P_Of_AA_Vs_Mu_S",sts,"AA",aas,"n",n,"k",k,"_phi",phi,"_",Codon.space,"_",initialization,".png")
image.title <- bquote("P("~alpha~") "~.(paste(sep="","Over ",maxT+1))~"Trajectories")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",naars,", T = ",ntrna,", S = ",sts,", n = ",n,", k = ",k,",")
image.subtitles[2] <- paste(sep="","= ",phi,", ",cs,", ",rate,", Initialized ",init,",")

png(image.name,height=1350,width=2200)
print(ggplot(q,aes(x=Fixation,y=q2))+geom_line(aes(color="black"))+geom_ribbon(data=q,aes(ymin=q1,ymax=q3),color="green",alpha=0.3)+geom_ribbon(data=q,aes(ymin=q0,ymax=q4),color="purple",alpha=0.2)+geom_ribbon(data=q,aes(ymin=qmin,ymax=qmax),color="red",alpha=0.1)+facet_grid(cols=vars(Mu),rows=vars(Amino_Acid))+labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+scale_color_manual("",values=c("Median"="black","Interquartile"="green","Inter-0.9-quantile"="purple","Extrema"="red"),labels=c("Median","Interquartile","Inter-0.9-quantile","Extrema"))+theme(legend.position="bottom",text=element_text(size=30),strip.text=element_text(size=15),strip.text.y=element_text(size=25),strip.text.x=element_text(size=25),axis.text.x=element_text(size=20))+ylab(bquote("P("~alpha~")"))+scale_x_continuous(breaks=seq(floor((maxfix+1)/4),3*floor((maxfix+1)/4),floor((maxfix+1)/4))))
dev.off()

#+geom_hline(yintercept=0.15)+geom_hline(yintercept=0.125)+geom_hline(yintercept=0.05)+geom_hline(yintercept=0.075)+geom_hline(yintercept=0.1)+geom_hline(yintercept=0.175)
