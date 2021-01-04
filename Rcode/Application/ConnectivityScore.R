rm(list=ls())
gc()

pparray = (1:12)[-c(5,6)]

pnodeattr = numeric(0)
pedgeattr = numeric(0)

for (pp in pparray){
    print(pp)
    load("../data/StringV10_ppi_pathway.Rdata")
    padat = ppi$pathwaydat
    padat[,2] = gsub("-","",padat[,2])
    pname = unique(padat[,1])[pp]
    proteins = toupper(padat[padat[,1] ==pname,2])
    genes =unlist(strsplit(padat[padat[,1] ==pname,3],split=", "))
    
    load("../data/TCGA_LUSC_CNA_Methylation_mRNA_RPPA.rda")
    sum(!colnames(RPPA)%in% toupper(padat[,2]))
    RPPA = RPPA[,colnames(RPPA)%in%proteins]
    mRNA = mRNA[,colnames(mRNA)%in%genes]
    CNA = CNA[,colnames(CNA)%in%genes]
    Methylation = Methylation[,colnames(Methylation)%in%genes]
    
    Y = cbind(Methylation,CNA,mRNA,RPPA)
    colnames(Y) = c(paste("Methyl_",colnames(Methylation),sep=""),paste("CNA_",colnames(CNA),sep=""),paste("mRNA_",colnames(mRNA),sep="")
    ,paste("RPPA_",colnames(RPPA),sep=""))
    
    
    p1 = ncol(Methylation)
    p2 = ncol(CNA)
    p3 = ncol(mRNA)
    p4 = ncol(RPPA)
    addr = cumsum(c(p1,p2,p3,p4))
    
    q=4
    
    chlist = vector("list",q)
    chlist[[1]] = 1:addr[1]
    chlist[[2]] = (addr[1]+1):addr[2]
    chlist[[3]] = (addr[2]+1):addr[3]
    chlist[[4]] = (addr[3]+1):addr[4]
    
    palist = vector("list",q)
    palist[[3]] = 1:addr[2]
    palist[[4]] = 1:addr[3]
    
    load(paste("../appl_results/LUSC_fit_pp",pp,".rda",sep=""))
    p = addr[q]
    
    B = matrix(0,p,p)
    cB = numeric(0)
    for (j in 1:q) {
        tmp = apply(fit[[j]]$eta[,,10001:30000],c(1,2),sum)/20000
        B[chlist[[j]],chlist[[j]]] = tmp
        cB = c(cB,tmp[lower.tri(tmp)])
        if (j%in%c(3,4)) {
            tmp = c(apply(fit[[j]]$Gamma[,,10001:30000],c(1,2),sum)/20000)
            B[as.matrix(expand.grid(chlist[[j]],palist[[j]]))] = tmp
            cB = c(cB,tmp)
        }
    }
    o = order(cB,decreasing=T)
    ii = max(which(cumsum(1-cB[o]) * (1/1:length(cB))<0.1))  ## FDR at 0.1
    cut = cB[o][ii]
    
    ########
    nodes = colnames(Y)
    matB = matrix(0,p,p)
    matB[B>=cut] = 1
    und.w = which(matB + t(matB)==2,arr.ind=T)
    und.w = unique(t(apply(und.w,1,sort)))
    dir.w = which(matB==1 & t(matB)==0,arr.ind=T)
    
    load(paste("../appl_results/LUSC_fit_str_pp",pp,".rda",sep=""))
    Coef = matrix(0,p,p)
    for (j in 1:q) {
        Coef[chlist[[j]],chlist[[j]]]= apply(fit[[j]]$A[,,10001:30000],c(1,2),sum)/20000
        if (j%in%c(3,4)) {
            Coef[as.matrix(expand.grid(chlist[[j]],palist[[j]]))] = c(apply(fit[[j]]$B[,,10001:30000],c(1,2),sum)/20000)
        }
    }
    
    undattr = cbind(nodes[und.w[,1]],nodes[und.w[,2]],"und",round(B[und.w],digit=2),sign(Coef[und.w]),rep(pname,nrow(und.w)))
    dirattr = cbind(nodes[dir.w[,2]],nodes[dir.w[,1]],"dir",round(B[dir.w],digit=2),sign(Coef[dir.w]),rep(pname,nrow(dir.w)))
    edgeattr = rbind(undattr,dirattr)
    onodes=nodes
    nodes = gsub("CNA_","CNA:",nodes)
    nodes = gsub("mRNA_","mRNA:",nodes)
    nodes = gsub("RPPA_","RPPA:",nodes)
    nodes = gsub("Methyl_","Methyl:",nodes)
    nodeattr = cbind(onodes,matrix(unlist(strsplit(nodes,split=":")),byrow=T,ncol=2))
    
    pedgeattr = rbind(pedgeattr,edgeattr)
    pnodeattr = rbind(pnodeattr,nodeattr)
}
pnodeattr = unique(pnodeattr)

write.table(pedgeattr,"edgeattr_LUSC_pw.txt",col.names=F,row.names=F,sep="\t")
write.table(pnodeattr,"nodeattr_LUSC_pw.txt",col.names=F,row.names=F,sep="\t")

##################### Compute pathway and platform specific connectivity scores ######
rm(list=ls())
gc()
dn.array = c("LUAD","LUSC","COAD","READ","UCEC","OV","SKCM")
pparray = (1:12)[-c(5,6)]

CS.all = vector("list",10)
for (pi in 1:length(pparray)) {
    pp = pparray[pi]


CS = matrix(nrow=9,ncol=length(dn.array))
for (dd in 1:length(dn.array)) {
    dn = dn.array[dd]

    ## pathway
    load("../data/StringV10_ppi_pathway.Rdata")
    padat = ppi$pathwaydat
    padat[,2] = gsub("-","",padat[,2])
    pname = unique(padat[,1])[pp]
    proteins = toupper(padat[padat[,1] ==pname,2])
    genes =unlist(strsplit(padat[padat[,1] ==pname,3],split=", "))

    ## Data
    load(paste("../data/TCGA_",dn,"_CNA_Methylation_mRNA_RPPA.rda",sep=""))
    sum(!colnames(RPPA)%in% toupper(padat[,2]))
    RPPA = RPPA[,colnames(RPPA)%in%proteins]
    mRNA = mRNA[,colnames(mRNA)%in%genes]
    CNA = CNA[,colnames(CNA)%in%genes]
    Methylation = Methylation[,colnames(Methylation)%in%genes]

    Y = cbind(Methylation,CNA,mRNA,RPPA)
    colnames(Y) = c(paste("Methyl_",colnames(Methylation),sep=""),paste("CNA_",colnames(CNA),sep=""),paste("mRNA_",colnames(mRNA),sep="")
    ,paste("RPPA_",colnames(RPPA),sep=""))


    p1 = ncol(Methylation)
    p2 = ncol(CNA)
    p3 = ncol(mRNA)
    p4 = ncol(RPPA)

    ## Edgeset
    edgeattr = read.table(paste("edgeattr_",dn,"_pw.txt",sep=""),header=F, sep = "\t")

    ## Compute Connectivity Score
    pedgeattr = edgeattr[edgeattr[,6]==pname,]

    # CNA-CNA #
    Tedges = (p2 * (p2-1))/2
    Nedges = length(intersect(grep("CNA",pedgeattr[,1]),grep("CNA",pedgeattr[,2])))
    CS[1,dd] = Nedges/Tedges
    
    # Methylation - Methylation #
    Tedges = (p1 * (p1-1))/2
    Nedges = length(intersect(grep("Methyl",pedgeattr[,1]),grep("Methyl",pedgeattr[,2])))
    CS[2,dd] = Nedges/Tedges
    
    # mRNA - mRNA #
    Tedges = (p3 * (p3-1))/2
    Nedges = length(intersect(grep("mRNA",pedgeattr[,1]),grep("mRNA",pedgeattr[,2])))
    CS[3,dd] = Nedges/Tedges
    
    # protein-protein #
    Tedges = (p4 * (p4-1))/2
    Nedges = length(intersect(grep("RPPA",pedgeattr[,1]),grep("RPPA",pedgeattr[,2])))
    CS[4,dd] = Nedges/Tedges
    
    # CNA -> mRNA #
    Tedges = p2 * p3
    Nedges = length(intersect(grep("CNA",pedgeattr[,1]),grep("mRNA",pedgeattr[,2])))
    CS[5,dd] = Nedges/Tedges
    
    # CNA -> protein #
    Tedges = p2 * p4
    Nedges = length(intersect(grep("CNA",pedgeattr[,1]),grep("RPPA",pedgeattr[,2])))
    CS[6,dd] = Nedges/Tedges
 
    # Methyl -> mRNA #
    Tedges = p1 * p3
    Nedges = length(intersect(grep("Methyl",pedgeattr[,1]),grep("mRNA",pedgeattr[,2])))
    CS[7,dd] = Nedges/Tedges

    # Methyl -> Protein #
    Tedges = p1 * p4
    Nedges = length(intersect(grep("Methyl",pedgeattr[,1]),grep("RPPA",pedgeattr[,2])))
    CS[8,dd] = Nedges/Tedges
    
    # mRNA -> Protein #
    Tedges = p3 * p4
    Nedges = length(intersect(grep("mRNA",pedgeattr[,1]),grep("RPPA",pedgeattr[,2])))
    CS[9,dd] = Nedges/Tedges
}

CS.all[[pi]] = CS

}

CSmat = do.call(rbind,CS.all)


### draw Heatmap ###
library(pheatmap)
library(RColorBrewer)
colnames(CSmat) = dn.array
rownames(CSmat) = paste(1:90,rep(c("CNA-CNA","Methylation-Methylation","mRNA-mRNA","Protein-Protein","CNA->mRNA","CNA->Protein","Methylation->mRNA","Methylation->Protein","mRNA->Protein"),10))
ann = data.frame(Pathway=rep(unique(padat[,1])[-c(5,6)],rep(9,10)))
rownames(ann) = rownames(CSmat)

colfunc1 <- colorRampPalette(c("lightgrey","red","black"))

row.lab=rep(c("CNA-CNA","Methylation-Methylation","mRNA-mRNA","Protein-Protein","CNA->mRNA","CNA->Protein","Methylation->mRNA","Methylation->Protein","mRNA->Protein"),10)
ann.colors =brewer.pal(10,"Paired")
names(ann.colors) =unique(padat[,1])[-c(5,6)]
ann.colors= list(Pathway=ann.colors)
pdf("CS_heatmap.pdf",height=7,width=4)
pheatmap(CSmat,color=colfunc1(200),cluster_cols=F,cluster_rows=F,fontsize_col=7,fontsize_row=5,annotation_row=ann,annotation_colors=ann.colors,labels_row=row.lab)
dev.off()

pdf("CS_sd1.pdf",height=2,width=7)
par(mar=c(0.5,4.5,0,1))
plot(apply(CSmat,1,sd)[45:1],type="h",pch=16,cex=1.2,lwd=5,col=rep(brewer.pal(10,"Paired"),rep(9,10))[45:1],xaxt="n",ylim=c(0,0.2))
dev.off()

pdf("CS_sd2.pdf",height=2,width=7)
par(mar=c(0.5,4.5,0,1))
plot(apply(CSmat,1,sd)[90:46],type="h",pch=16,cex=1.2,lwd=5,col=rep(brewer.pal(10,"Paired"),rep(9,10))[90:46],xaxt="n",ylim=c(0,0.2))
dev.off()


