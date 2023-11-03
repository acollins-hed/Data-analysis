library("tidyverse")
library("ggsankey")
library("dplyr")

codon.frequencies <- c("all","frequencies are 10")
rate <- "rate dep"
initialization <- "all_0s"
init <- "all 0s"
codon.space <- "Cring"
cs <- "codon ring"
phi <- 0.5
k <- 4
n <- 8

MU <- c("0.0001","0.001","0.005","0.01","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2")

Mu <- MU[1]

file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_code.dat")

print(paste(sep="","Loading file ",file))
code <- read.table(file,header=TRUE,row.names=NULL)
print(paste(sep="",file," loaded"))

file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_codon.dat")

print(paste(sep="","Loading file ",file))
codon <- read.table(file,header=TRUE,row.names=NULL)
print(paste(sep="",file," loaded"))

sts <- unique(codon$Site_Type)
nst <- length(sts)
aarss <- unique(code$aaRS)
trnas <- unique(code$tRNA)
naars <- length(aarss)
ntrna <- length(trnas)
maxfix <- max(code$Fixation)
maxT <- max(code$Trajectory)

code.temp <- code[which(code$Trajectory == 0 & code$Fixation == maxfix),]
codon.temp <- codon[which(codon$Trajectory == 0 & codon$Fixation == maxfix),]

test <- cbind(codon.temp[,1:4],codon.temp[,4:5])
colnames(test)[5] <- "tRNA"

temp <- test[which(test$tRNA == trnas[1]),]
temp1 <- code.temp[which(code.temp$tRNA == trnas[1]),3:5]
temp2 <- temp[1,]
for(i in 2:dim(temp1)[1]){
    temp2 <- rbind(temp2,temp[1,])
}
temp3 <- cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
for(j in 2:dim(temp)[1]){
    temp2 <- temp[j,]
    for(i in 2:dim(temp1)[1]){
        temp2 <- rbind(temp2,temp[j,])
    }
    temp4 <- cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
    temp3 <- rbind(temp3,temp4)
}
forsankey <- temp3 # cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
colnames(forsankey)[8] <- "Codon_Frequency"

for(t in trnas[-1]){
    temp <- test[which(test$tRNA == t),]
    temp1 <- code.temp[which(code.temp$tRNA == t),3:5]
    temp2 <- temp[1,]
    for(i in 2:dim(temp1)[1]){
        temp2 <- rbind(temp2,temp[1,])
    }
    temp3 <- cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
    for(j in 2:dim(temp)[1]){
        temp2 <- temp[j,]
        for(i in 2:dim(temp1)[1]){
            temp2 <- rbind(temp2,temp[j,])
        }
        if(j == 2){
            temp5 <- temp2
        }else{
            temp5 <- rbind(temp5,temp2)
        }
    }
    temp4 <- cbind(temp5[,1:5],temp1[,2:3],temp5[,6])
    colnames(temp4)[8] <- "Codon_Frequency"
    colnames(temp3)[8] <- "Codon_Frequency"
    temp3 <- rbind(temp3,temp4)
    #temp2 <- cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
    colnames(temp3)[8] <- "Codon_Frequency"
    forsankey <- rbind(forsankey,temp3)
}

forsankey <- cbind(forsankey[,1:6],forsankey[,6:8])
forsankey <- cbind(forsankey,forsankey$Prob_Interaction*forsankey$Codon_Frequency)
colnames(forsankey)[c(7,10)] <- c("Amino_Acid","Weight")

df <- forsankey %>% make_long(Codon,tRNA,aaRS,Amino_Acid,Site_Type,value=Weight)

image.name <- paste(sep="","Sankey_Diagram_A",naars,"T",ntrna,"n",n,"k",k,"_phi",phi,"Mu",Mu,"_",codon.space,"_",initialization,".png")
image.title <- bquote("Sankey Diagram for Final Fixation "~.(paste(sep="","Afterer ",maxfix+1))~"Fixations")
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",naars,", T = ",ntrna,", S = ",nst,", n = ",n,", k = ",k,",")
image.subtitles[2] <- paste(sep="","= ",phi,", Mu = ",Mu,", ",cs,", ",rate,", Initialized ",init,",")

png(image.name,height=1350,width=2200)
print(ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node, value=value)) +geom_sankey(flow.alpha = .6,node.color = "gray30") +geom_sankey_label(size = 10, color = "white", fill = "gray40") +scale_fill_viridis_d() +theme_sankey(base_size = 24) +labs(x = NULL) +theme(legend.position = "none") +labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+geom_sankey_label(size = 20, color = 1, fill = "white")+theme(text = element_text(size=42),axis.text.x=element_text(size=50)))
dev.off()

df <- cbind(df,rep(Mu,dim(df)[1]))
colnames(df)[dim(df)[2]] <- "Mu"

for(Mu in MU[-1]){

    file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_code.dat")

    print(paste(sep="","Loading file ",file))
    code <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," loaded"))
    
    file <- paste(sep="","/media/sciship/Seagate Expansion Drive/data/A5T5n8k4_phi0.5_all0s_St_Freq_All_20/atINFLATE_A5T5S5n8k4phi0.5Mu",Mu,"fix1000C2all0_rate_dep_uniform_aas_codon.dat")
    
    print(paste(sep="","Loading file ",file))
    codon <- read.table(file,header=TRUE,row.names=NULL)
    print(paste(sep="",file," loaded"))
    
    sts <- unique(codon$Site_Type)
    nst <- length(sts)
    aarss <- unique(code$aaRS)
    trnas <- unique(code$tRNA)
    naars <- length(aarss)
    ntrna <- length(trnas)
    maxfix <- max(code$Fixation)
    maxT <- max(code$Trajectory)
    
    code.temp <- code[which(code$Trajectory == 0 & code$Fixation == maxfix),]
    codon.temp <- codon[which(codon$Trajectory == 0 & codon$Fixation == maxfix),]
    
    test <- cbind(codon.temp[,1:4],codon.temp[,4:5])
    colnames(test)[5] <- "tRNA"
    
    temp <- test[which(test$tRNA == trnas[1]),]
    temp1 <- code.temp[which(code.temp$tRNA == trnas[1]),3:5]
    temp2 <- temp[1,]
    for(i in 2:dim(temp1)[1]){
        temp2 <- rbind(temp2,temp[1,])
    }
    temp3 <- cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
    for(j in 2:dim(temp)[1]){
        temp2 <- temp[j,]
        for(i in 2:dim(temp1)[1]){
            temp2 <- rbind(temp2,temp[j,])
        }
        temp4 <- cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
        temp3 <- rbind(temp3,temp4)
    }
    forsankey <- temp3 # cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
    colnames(forsankey)[8] <- "Codon_Frequency"
    
    for(t in trnas[-1]){
        temp <- test[which(test$tRNA == t),]
        temp1 <- code.temp[which(code.temp$tRNA == t),3:5]
        temp2 <- temp[1,]
        for(i in 2:dim(temp1)[1]){
            temp2 <- rbind(temp2,temp[1,])
        }
        temp3 <- cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
        for(j in 2:dim(temp)[1]){
            temp2 <- temp[j,]
            for(i in 2:dim(temp1)[1]){
                temp2 <- rbind(temp2,temp[j,])
            }
            if(j == 2){
                temp5 <- temp2
            }else{
                temp5 <- rbind(temp5,temp2)
            }
        }
        temp4 <- cbind(temp5[,1:5],temp1[,2:3],temp5[,6])
        colnames(temp4)[8] <- "Codon_Frequency"
        colnames(temp3)[8] <- "Codon_Frequency"
        temp3 <- rbind(temp3,temp4)
                                        #temp2 <- cbind(temp2[,1:5],temp1[,2:3],temp2[,6])
        colnames(temp3)[8] <- "Codon_Frequency"
        forsankey <- rbind(forsankey,temp3)
    }
    
    forsankey <- cbind(forsankey[,1:6],forsankey[,6:8])
    forsankey <- cbind(forsankey,forsankey$Prob_Interaction*forsankey$Codon_Frequency)
    colnames(forsankey)[c(7,10)] <- c("Amino_Acid","Weight")
    
    temp.df <- forsankey %>% make_long(Codon,tRNA,aaRS,Amino_Acid,Site_Type,value=Weight)
    
    image.name <- paste(sep="","Sankey_Diagram_A",naars,"T",ntrna,"n",n,"k",k,"_phi",phi,"Mu",Mu,"_",codon.space,"_",initialization,".png")
    image.title <- bquote("Sankey Diagram for Final Fixation "~.(paste(sep="","Over ",maxT+1))~"Trajectories,"~.(paste(sep="",maxfix+1," Fixations Each")))
    image.subtitles <- c(0)
    image.subtitles[1] <- paste(sep="","A = ",naars,", T = ",ntrna,", S = ",nst,", n = ",n,", k = ",k,",")
    image.subtitles[2] <- paste(sep="","= ",phi,", Mu = ",Mu,", ",cs,", ",rate,", Initialized ",init,",")
    
    png(image.name,height=1350,width=2200)
    print(ggplot(temp.df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node, value=value)) +geom_sankey(flow.alpha = .6,node.color = "gray30") +geom_sankey_label(size = 10, color = "white", fill = "gray40") +scale_fill_viridis_d() +theme_sankey(base_size = 24) +labs(x = NULL) +theme(legend.position = "none") +labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+geom_sankey_label(size = 20, color = 1, fill = "white")+theme(text = element_text(size=42),axis.text.x=element_text(size=50)))
    dev.off()

    temp.df <- cbind(temp.df,rep(Mu,dim(temp.df)[1]))
    colnames(temp.df)[dim(temp.df)[2]] <- "Mu"
    df <- rbind(df,temp.df)
}

df$Mu <- paste(sep="","Mu = ",df$Mu)

image.name <- paste(sep="","Sankey_Diagram_multipanel_A",naars,"T",ntrna,"n",n,"k",k,"_phi",phi,"_",codon.space,"_",initialization,".png")
image.title <- bquote("Sankey Diagram for Final Fixation "~.(paste(sep="","Over ",maxT+1))~"Trajectories,"~.(paste(sep="",maxfix+1," Fixations Each")))
image.subtitles <- c(0)
image.subtitles[1] <- paste(sep="","A = ",naars,", T = ",ntrna,", S = ",nst,", n = ",n,", k = ",k,",")
image.subtitles[2] <- paste(sep="","= ",phi,", ",cs,", ",rate,", Initialized ",init,",")

png(image.name,height=1350,width=2200)
print(ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node, value=value)) +geom_sankey(flow.alpha = .6,node.color = "gray30") +geom_sankey_label(size = 8, color = "white", fill = "gray40") +scale_fill_viridis_d() +theme_sankey(base_size = 20) +labs(x = NULL) +theme(legend.position = "none") +labs(title=image.title,subtitle=bquote(.(image.subtitles[1])~phi~.(image.subtitles[2])~.(codon.frequencies[1])~s[alpha]~.(codon.frequencies[2])))+geom_sankey_label(size = 8, color = 1, fill = "white")+theme(text = element_text(size=46),axis.text.x=element_text(size=20))+facet_wrap(~Mu))
dev.off()

write.table(df,file="A5T5S5_phi0.5_all0_rate_dep_cring_sankey.dat",row.names=FALSE)
