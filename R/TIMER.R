# using TIMER core algo to NAR server 
# data.bulk is the bulk data.
# bug report: wnchang@iu.edu
TIMER_bypassTCGA <- function(data.matrix = data.bulk){
	#cc = tolower(cc)
	#cancer_lib <- c("BLCA", "BRCA", "CESC", "COAD", "DLBC","ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", 
	#			"LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "UCS", "UVM")
	#cancer_lib_lll <- tolower(cancer_lib)
	#if(!cc%in%cancer_lib_lll)
	#	stop("please input the cancer type which exist in TCGA library. For example BLCA, COAD, SKCM...")
	#if(cc=='skcm')cc.type='06A' else cc.type='01A'

	##----- setup parameters and establish the output file -----##
	signature.genes=c('CD19','TRAT1','CD8B','CCR3','CD163','CCL17')
	names(signature.genes)=c('B_cell','T_cell.CD4','T_cell.CD8','Neutrophil','Macrophage','DC')

	##----- load and process gene expression data -----##
	#check sample name
	if(length(colnames(data.matrix)) == 0) {
	warning("input data do NOT have colnames")
	colnames(data.matrix) <- paste( "Setsample", 1:ncol(data.matrix), sep="")    
	}
	data.matrix <- rm_zero_row(data.matrix)
	#data_ttt <- filter_gene_name(data_t)
	data_ttt <- data.matrix
	dd <- data_ttt
	#dd <- get(load(file_str) )
	dd=as.matrix(dd)
	mode(dd)='numeric'	
	#if(!cc %in% c('gbm','ov','esca','stad'))dd=dd*1e6   ## rsem scaled estimates needs multiply 1e6, Array or RPKM does not need.
	tmp=strsplit(rownames(dd),'\\|')
	tmp=sapply(tmp,function(x)x[[1]])
	tmp.vv=which(nchar(tmp)>1)
	rownames(dd)=tmp
	dd=dd[tmp.vv,]

	#chang: load TIMER.Rdata from local directory for replacing all load opearation below.
	#load("./TIMER.RData")
	##----- load immune marker genes from Abbas et al., 2005 -----##
	#tmp=read.csv("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/IRIS-marker-gene.txt", ,header=T,sep='\t',stringsAsFactors=F) # * add / before data
	#IRIS_marker_gene = tmp
	#save(IRIS_marker_gene,file="IRIS_marker_gene.RData")
	tmp = IRIS_marker_gene		#chang : save csv into RData
	marker.list=tmp[,1]
	names(marker.list)=tmp[,7]
	names(marker.list)=gsub(' ','_',tmp[,7])
	names(marker.list)=gsub('Dendritic_Cell','DC',names(marker.list))
	names(marker.list)=gsub('Neutrophil','Neutrophils',names(marker.list))
	gsub('Cell','cell',names(marker.list))->names(marker.list)

	##----- load reference data of sorted immune cells -----##
	#load("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/HPCTimmune.Rdata") # * add / before data， load HPCT.immune.

	##----- load and process tumor purity data -----##
	#AGP=read.table(paste('C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/AGP/AGP-',cc,'.txt',sep=''),sep='\t',header=T) # * add / before datas
	#AGP=AGP[which(AGP[,'PoP']>0.01),]
	tmp=strsplit(rownames(HPCT.immune),';')
	AffyIDtoGenes=sapply(tmp,function(x)x[[1]])
	names(AffyIDtoGenes)=sapply(tmp,function(x)x[[2]])
	marker.list.genes=AffyIDtoGenes[marker.list]

	##----- function to edit TCGA ID, with the option of keeping the first num.res fields -----##
	#getID <- function(sID,num.res=3){
 	#	mm=c()
  	#	for(id in sID){
   	#		tmp=unlist(strsplit(id,'-'))
    #		if(length(tmp)==1){
    #	 		tmp=unlist(strsplit(id,'\\.'))
    #		}
    #		ll='TCGA'
    #		for(j in 2:num.res){
    #			ll=paste(ll,tmp[j],sep='-')
   	#	 	}
    #		mm=c(mm,ll)
 	#	}
  	#	return(mm)
	#}
	#rownames(AGP)=getID(AGP[,1],4)
	#colnames(dd)=getID(colnames(dd),4)


	##----- Select single reference samples of pre-selected immune cell types -----##
	B_cell=362:385
	T_cell.CD4=grep('T_cell.CD4',colnames(HPCT.immune))
	T_cell.CD8=grep('T_cell.CD8',colnames(HPCT.immune))
	NK=328:331
	Neutrophil=344:361
	Macrophage=66:80
	DC=151:238
	curated.ref=HPCT.immune[,c(B_cell,T_cell.CD4,T_cell.CD8,NK,Neutrophil,Macrophage,DC)]

	curated.cell.types=colnames(curated.ref)
	names(curated.cell.types)=c(rep('B_cell',length(B_cell)),rep('T_cell.CD4',length(T_cell.CD4)),rep('T_cell.CD8',length(T_cell.CD8)),rep('NK',length(NK)),rep('Neutrophil',length(Neutrophil)),rep('Macrophage',length(Macrophage)),rep('DC',length(DC)))

	#load('C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/curated.ref.genes.Rdata')	#chang: curated.ref.genes 

	##----- Combine TCGA gene expression profiles with the selected reference data, remove batch effect and aggregate samples of each immune category by taking the median -----##
	RemoveBatchEffect <- function(){
  		library(sva)
  		tmp.dd=as.matrix(dd)
  		tmp=sapply(strsplit(rownames(dd),'\\|'),function(x)x[[1]])
  		rownames(tmp.dd)=tmp
  		tmp.dd=tmp.dd[which(nchar(tmp)>1),]
  		tmp.ss=intersect(rownames(tmp.dd),rownames(curated.ref.genes))
  		N1=ncol(tmp.dd)
  		tmp.dd=cbind(tmp.dd[tmp.ss,],curated.ref.genes[tmp.ss,])
  		tmp.dd=as.matrix(tmp.dd)
  		mode(tmp.dd)='numeric'
  		N2=ncol(curated.ref.genes)
  		tmp.batch=c(rep(1,N1),rep(2,N2))
  		tmp.dd0=ComBat(tmp.dd,tmp.batch,c())
  		dd.br=tmp.dd0[,1:N1]
  		curated.ref.genes.br=tmp.dd0[,(N1+1):(N1+N2)]
  		tmp0=c()
  		for(kk in unique(names(curated.cell.types))){
   		 tmp.vv=which(names(curated.cell.types)==kk)
   		 tmp0=cbind(tmp0,apply(curated.ref.genes.br[,tmp.vv],1,median,na.rm=T))
  		}
  		curated.ref.genes.agg.br=tmp0
  		colnames(curated.ref.genes.agg.br)=unique(names(curated.cell.types))
  		#rownames(curated.ref.genes.agg.br)=rownames(curated.ref.genes.br)
  		return(list(dd=dd.br,rr=curated.ref.genes.br,rrg=curated.ref.genes.agg.br))
	}

	tmp=RemoveBatchEffect()
	dd.br=tmp$dd
	curated.ref.genes.br=tmp$rr
	curated.ref.genes.agg.br=tmp$rrg


	##----- function to calculate the residuals from regression -----##
	fn <- function(beta0,XX,Y)return(log(sum(abs(Y-XX%*%beta0))))

	##----- function to select genes with expression values negatively correlated with tumor purity -----##
	getPurityGenes <- function(dd){
	#getPurityGenes <- function(dd,AGP,thr.p=0.05,thr.c=0,mode='env'){
		#	 tmp.ss=intersect(colnames(dd),rownames(AGP))
		# 	if(length(tmp.ss)==0){
		#   	colnames(dd)=getID(colnames(dd))
		#  		tmp.ss=intersect(colnames(dd),rownames(AGP))
		# 	}
		# 	tmp.dd=dd[,tmp.ss]	#original
		# 	#tmp.dd <- dd    #w
		# 	tmp=lapply(rownames(tmp.dd),function(x)cor.test(tmp.dd[x,],as.numeric(AGP[colnames(tmp.dd),2]),method='s'))
		# 	tmp.pp=sapply(tmp,function(x)x$p.value)
		# 	tmp.cor=sapply(tmp,function(x)x$estimate)
		# 	names(tmp.pp)=names(tmp.cor)=rownames(dd)
		# 	if(mode=='env')vv=names(which(tmp.pp <=thr.p&tmp.cor < thr.c))
		# 	if(mode=='tumor')vv=names(which(tmp.pp <=thr.p&tmp.cor > thr.c))
		#	return(vv)

 		#single cell data do NOT have overlap sample with AGP, thus we assume all genes are purity genes.(chang.)
 		return(rownames(dd))
	}

	##----- selection genes negatively correlated with purity and overlap with immune marker genes -----##
	#vv.t=getPurityGenes(dd,AGP,thr.p=0.05,thr.c= -0.2)
	vv.t = getPurityGenes(dd)
	vv.t=intersect(vv.t,rownames(curated.ref.genes.agg.br))
	vv=intersect(vv.t,marker.list.genes)

	##----- remove outlier genes whose expression may drive the colinearity of similar covariates in the regression -----##
	RemoveOutliers <- function(vv, ref.dd, thr.q=0.99){
 	 ## removes upper thr.q quantile for every reference feature
 	 remove.vv=c()
 	 for(i in 1:ncol(ref.dd)){
  	  tmp=quantile(ref.dd[vv,i],thr.q)[1]
  	  tmp.vv=which(ref.dd[vv,i]>tmp)
 	  remove.vv=c(remove.vv,tmp.vv)
	 }
 	 remove.vv=unique(remove.vv)
 	 return(vv[-remove.vv])
	}

	##---- calculate differences between the correlations of reference immune cells using Pearson's or Spearman's correlations -----##
	tmp.diff=sum(sum(abs(cor(curated.ref.genes.agg.br[vv,],method='p')-cor(curated.ref.genes.agg.br[vv,],method='s'))))

	if(tmp.diff>= -10000){
 		vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4])
 		 vv=vv0
	}

	#cat("Number of genes inversely correlated with purity is ",length(vv.t),'\n\n',sep='',file='output-statistics.txt',append=T)
	#cat("Number of immune genes inversely correlated with purity is ",length(vv),'\n\n',sep='',file='output-statistics.txt',append=T)

	##----- calculate the significance of enrichment for purity selected genes to immune marker genes -----##
	tmp.ss0=intersect(rownames(curated.ref.genes.agg.br),rownames(dd.br))
	n.immune=length(intersect(marker.list.genes,tmp.ss0))
	cat("Test if immune genes are enriched for inverse correlation with purity: \n\n",file='output-statistics.txt',append=T)
	sink(file='output-statistics.txt',append=T);print(fisher.test(matrix(c(length(vv),length(vv.t)-length(vv),n.immune,length(tmp.ss0)-n.immune),2,2)));sink()


	##----- function to process deconvolution method in batch -----##
	BatchFractions <- function(XX,YYd){
  		Fmat=c()
  		for(i in 1:ncol(YYd)){
  			YY=YYd[,i]
   			tmp.F=getFractions.Abbas(XX,YY)
    		#tmp.F=getFractions.Optim(XX,YY)
    		Fmat=rbind(Fmat,tmp.F)
  		}
  		rownames(Fmat)=colnames(YYd)
  		colnames(Fmat)=colnames(XX)
  		return(Fmat)
	}

	##----- perform batch deconvolution -----##
	XX=curated.ref.genes.agg.br[vv,c(-4)]	#chang note: delete NK cell during deconvolution
	YYd=dd.br[vv,]
	Fmat=BatchFractions(XX,YYd)

	if( !is.na( cor(Fmat[,2],Fmat[,3])[1] ) ){
		print("Remove co-linear!!!")
		##----- CD4 and CD8 T cells are likely to be similar, resulting in colinearity. Codes below are procedures to remove outlier genes that may result in colinearity until the two covariates are linearly separable. -----##
		if(cor(Fmat[,2],Fmat[,3])<= -0.2){
			if(tmp.diff>=1){
				tmp.cor=c()
				thr.qlist=c(0.99)
				for(tq in thr.qlist){
				 vv=intersect(vv.t,marker.list.genes)
				 vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4],tq)
				 vv=vv0
				 XX=curated.ref.genes.agg.br[vv,]
				 YYd=dd.br[vv,]
				 tmp.Fmat=BatchFractions(XX,YYd)
				 tmp.cor=c(tmp.cor,cor(tmp.Fmat[,2],tmp.Fmat[,3],method='s'))
				}
				tmp.vv=which.max(tmp.cor)
				tq=thr.qlist[tmp.vv]
				vv=intersect(vv.t,marker.list.genes)
				vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,c(-4)],tq)
				vv=vv0
				XX=curated.ref.genes.agg.br[vv,c(-4)]
				YYd=dd.br[vv,]
				Fmat=BatchFractions(XX,YYd)
				#Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]	#chang: delete cc.type
				#rownames(Fmat0.p)=getID(rownames(Fmat0.p))
			}
		}

		count=0
		while(cor(Fmat[,2],Fmat[,3])<=-0.3){
			print("Remove co-linearality iterativelly!!~")
			count=count+1
			#print(count)
		 if(length(vv)<=50)break
		 vv=vv[-as.numeric(names(table(apply(dd[vv,],2,which.max))))]
		 XX=curated.ref.genes.agg.br[vv,c(-4)]
		 #XX=XX[,!colnames(XX) %in% tmp.remove]
		 YYd=dd.br[vv,]
		 Fmat=BatchFractions(XX,YYd)
			if(count>=30)break
		}
	}
	#Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]
	#rownames(Fmat0.p)=getID(rownames(Fmat0.p))

	TIMER_fraction <- Fmat
	X <- YYd	# chang: bulk data subset
	S <- XX #negative value	# signature matrix
	return(list(Timer_fraction = TIMER_fraction, X = X, S = S))

}

##----- Constrained regression method implemented in Abbas et al., 2009 -----##
getFractions.Abbas <- function(XX, YY, w=NA){
  ## XX is immune expression data
  ## YY is cancer expression data
  ss.remove=c()
  ss.names=colnames(XX)
  while(T){
    if(length(ss.remove)==0)tmp.XX=XX else{
      if(is.null(ncol(tmp.XX)))return(rep(0, ncol(XX)))
      tmp.XX=tmp.XX[, -ss.remove]
    }
    if(length(ss.remove)>0){
      ss.names=ss.names[-ss.remove]
      if(length(ss.names)==0)return(rep(0, ncol(XX)))
    }
    if(is.na(w[1]))tmp=lsfit(tmp.XX, YY, intercept=F) else tmp=lsfit(tmp.XX, YY, w, intercept=F)
    if(is.null(ncol(tmp.XX)))tmp.beta=tmp$coefficients[1] else tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
    if(min(tmp.beta>0))break
    ss.remove=which.min(tmp.beta)
  }
  tmp.F=rep(0, ncol(XX))
  names(tmp.F)=colnames(XX)
  tmp.F[ss.names]=tmp.beta
  return(tmp.F)
}



extract_TIMER_sigMat <- function()
{
	tmp=strsplit(rownames(HPCT.immune),';')
    AffyIDtoGenes=sapply(tmp,function(x)x[[1]])
    rownames(HPCT.immune) <- AffyIDtoGenes
	
	B_cell=362:385
	T_cell.CD4=grep('T_cell.CD4',colnames(HPCT.immune))
	T_cell.CD8=grep('T_cell.CD8',colnames(HPCT.immune))
	NK=328:331
	Neutrophil=344:361
	Macrophage=66:80
	DC=151:238
	curated.ref=HPCT.immune[,c(B_cell,T_cell.CD4,T_cell.CD8,NK,Neutrophil,Macrophage,DC)]

	curated.cell.types=colnames(curated.ref)
	names(curated.cell.types)=c(rep('B_cell',length(B_cell)),rep('T_cell.CD4',length(T_cell.CD4)),rep('T_cell.CD8',length(T_cell.CD8)),rep('NK',length(NK)),rep('Neutrophil',length(Neutrophil)),rep('Macrophage',length(Macrophage)),rep('DC',length(DC)))
	
	tmp0=c()
  	for(kk in unique(names(curated.cell.types))){
   		tmp.vv=which(names(curated.cell.types)==kk)
   		tmp0=cbind(tmp0,apply(curated.ref.genes[,tmp.vv],1,median,na.rm=T))			#chang: change to curated.ref.genes
  	}
  	curated.ref.genes.agg.br=tmp0
  	colnames(curated.ref.genes.agg.br)=unique(names(curated.cell.types))

	XX=curated.ref.genes.agg.br[,c(-4)]	#remove NK cell

	return(XX)
}
