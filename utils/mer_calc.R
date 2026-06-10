# Read in arguments
args <- commandArgs(TRUE) 
fl <- args[1]
pedigree <- args[2]
chr <- args[3]
tmpdir <- args[4]
outdir <- args[5]

cat(
   "fl -> ", fl, "\n",
   "pedigree -> ", pedigree, "\n",
   "chr -> ", chr, "\n",
   "tmpdir -> ", tmpdir, "\n",
   "outdir -> ", outdir, "\n"
)


#######################
## mdl2f Function
#######################

mdl2f<- function(geno, genoF0){
   genoF0<- genoF0[, -c(1:2)]
      genoF0<- as.matrix(genoF0)
   ids<- colnames(geno)[-c(1:2)]
   geno<- geno[, -c(1:2)]
      geno<- as.matrix(geno)
   idx<-    sweep(geno!=0, 1, genoF0[,1]==0 & genoF0[,2]==0, "&")
      idx<- sweep(geno==2, 1, genoF0[,1]==0 & genoF0[,2]==1, "&") | idx
      idx<- sweep(geno!=1, 1, genoF0[,1]==0 & genoF0[,2]==2, "&") | idx
      idx<- sweep(geno==2, 1, genoF0[,1]==1 & genoF0[,2]==0, "&") | idx
      idx<- sweep(geno==0, 1, genoF0[,1]==1 & genoF0[,2]==2, "&") | idx
      idx<- sweep(geno!=1, 1, genoF0[,1]==2 & genoF0[,2]==0, "&") | idx
      idx<- sweep(geno==0, 1, genoF0[,1]==2 & genoF0[,2]==1, "&") | idx
      idx<- sweep(geno!=2, 1, genoF0[,1]==2 & genoF0[,2]==2, "&") | idx

   idx
}


#######################
## mdl2f.m Function
#######################

mdl2f.m<- function(geno, genoF0){
   genoF0<- genoF0[, -c(1:2)] # (father, mother)
      genoF0<- as.matrix(genoF0)
   ids<- colnames(geno)[-c(1:2)]
   geno<- geno[, -c(1:2)]
      geno<- as.matrix(geno)
   idx<-    geno==1
      idx<- sweep(geno!=0, 1, genoF0[,2]==0, "&") | idx
      idx<- sweep(geno!=2, 1, genoF0[,2]==2, "&") | idx

   idx
}


#######################
## mdlf Function
## Calls on mdl2f and mdl2f.m Functions

# Input Parameters:   
# fl: vcf file
# ped: pedigree (id, father/sire, mother/dam, ...)
# tmpdir: directory to save temporary data
# outdir: directory to save chromosome-specific csv files

# Output Files: 
# imr_chr.csv - per-sample F_MISS
# lmr_chr.csv - per-SNP F_MISS
# ier_chr.csv - per-sample MER
# ler_chr.csv - per-SNP MER
#######################


mdlf<- function(fl, ped, outdir, chr, tmpdir="Tmp"){
   tdC<- FALSE
   if(!dir.exists(tmpdir)){
      dir.create(tmpdir, recursive=TRUE)
      tdC<- TRUE
      cat("Folder \"", tmpdir, "\" created\n", sep="")
   }
   tg<- sapply(strsplit(fl, "/"), function(x) x[length(x)])
      tg<- sapply(strsplit(tg, "\\."), function(x) x[1])
   l<- system(paste("zcat ", fl, " | grep -n -m 1 CHROM | cut -d \":\" -f 1", sep=""), intern=TRUE)
      l<- as.integer(l)

   hd<- read.table(fl, nrows=1, skip=l-1, header=FALSE, comment.char="@", colClasses = "character")

   idsT<- t(hd[,-(1:9)])[,1]
      names(idsT)<- NULL
      length(idsT)
   # must in pedigree
   cns<- colnames(ped)
   if(!is.element("id", tolower(cns))){
      idx<- match("sample_id", tolower(cns)) # formerly match("rfid")
      if(is.na(idx)){
         stop("ped: ID column name should be given as 'id' or 'rfid'", call.=FALSE)
      }else colnames(ped)[idx]<- "id"
   }else{
      idx<- match("id", tolower(cns))
      colnames(ped)[idx]<- "id"
   }
   if(!is.element("father", tolower(cns))){
      idx<- match("sire_sample_id", tolower(cns))
      if(is.na(idx)){
         stop("ped: paternal column name should be given as 'father' or 'sire'", call.=FALSE)
      }else colnames(ped)[idx]<- "father"
   }else{
      idx<- match("father", tolower(cns))
      colnames(ped)[idx]<- "father"
   }
   if(!is.element("mother", tolower(cns))){
      idx<- match("dam_sample_id", tolower(cns))
      if(is.na(idx)){
         stop("ped: maternal column name should be given as 'mother' or 'dam'", call.=FALSE)
      }else colnames(ped)[idx]<- "mother"
   }else{
      idx<- match("mother", tolower(cns))
      colnames(ped)[idx]<- "mother"
   }
   idx<- is.element(ped$mother,idsT) | is.element(ped$father,idsT) 
   ids<- unique(c(ped$id[idx], ped$mother[idx], ped$father[idx])) 
      ids<- intersect(ids, idsT)
      length(ids)

   idx<- match(ids, ped$id)
   idxF<- is.element(ped$father[idx], idsT)
   idxM<- is.element(ped$mother[idx], idsT)
   idsF<- ids[idxF & !idxM]
   idsM<- ids[!idxF & idxM]
   idsB<- ids[idxF & idxM]
   fms<- paste(ped$father, ped$mother, sep="_")
      fms<- fms[match(idsB, ped$id)]
      fms<- unique(fms)
   cat("[",length(idsF), ",", length(idsM), ",", length(idsB), "] samples have father, mother or both parents\n", sep="")
   cat("Those having both parents are from ",length(fms), " families\n", sep="")
   rm(idx, idxF, idxM, idsF, idsM, fms, l, hd, idsT)
   
   of<- file.path(tmpdir, paste(tg, "_gt_", chr, ".vcf", sep=""))

   str<- paste(ids, collapse=",")
      str<- paste(bcftools, " query -s ", str, " -r ", chr, sep="")
      str<- paste(str, " -f '%CHROM\t%POS\t[ %GT]\n' -o ", of, " ", fl, sep="")
   cat("\n")
   cat(str)
   cat("\n")
   system(str)
  
   gt<- read.table(of, header=FALSE, na.strings=c(".", "./."))
      rownames(gt)<- paste(gt[,1], gt[,2], sep="_")
      colnames(gt)<- c("Chr", "Pos", ids)

   geno<- gt[,-c(1:2)]
      geno<- (geno == "0/1")*1 + (geno == "1/1")*2
      geno<- cbind(gt[,c(1:2)], geno)

   rm(gt)
   system(paste("rm ", of, sep=""))

   nSnp<- nrow(geno)
   
   idx<- match(idsB, ped$id)
   pF<- ped$father[idx]
   pM<- ped$mother[idx]
   pB<- cbind(pF, pM, paste(pF, pM, sep="_"))
      pB<- pB[!duplicated(pB[,3]),]
   mdlB<- matrix(NA, nrow=nSnp, ncol=length(idsB))
      rownames(mdlB)<- paste(geno[,1], geno[,2], sep="_")
      colnames(mdlB)<- idsB

   print(pB)
   print(dim(pB))

   for(n in 1:nrow(pB)){
      idsT<- ped$id[ped$father == pB[n,1] & ped$mother == pB[n,2]]
         idsT<- intersect(idsT, idsB)
      gtT<- geno[, c(1:2, match(idsT,colnames(geno)))]
      gtP<- geno[, c(1:2, match(pB[n,1:2],colnames(geno)))]
      jj<- match(idsT, colnames(mdlB))
      if(chr == "chrX"){
         idx<- match(idsT, ped$id)
         sx<- toupper(ped$sex[idx])
         idx<- sx=="M" | sx=="MALE"
         if(any(idx)){
            jt<- !is.element(colnames(gtT), idsT[!idx])
            mdlB[,jj[idx]]<- mdl2f.m(geno=gtT[,jt], genoF0=gtP)
         }
         if(any(!idx)){
            jt<- !is.element(colnames(gtT), idsT[idx])
            mdlB[,jj[!idx]]<- mdl2f(geno=gtT[,jt], genoF0=gtP)
         }
         rm(idx)
      }else{
         mdlB[,jj]<- mdl2f(geno=gtT, genoF0=gtP)
      }

      rm(n, idsT, gtT, gtP)
   }

   imr<- apply(is.na(mdlB), 2, mean)
      imr<- data.frame(id=names(imr), r=imr)
      rownames(imr)<- NULL
   lmr<- apply(is.na(mdlB), 1, mean)
      st<- strsplit(names(lmr), "_")
      st<- data.frame(chr=sapply(st, function(x)x[1]), pos=sapply(st,function(x)as.integer(x[2])))
      lmr<- cbind(st, r=lmr)
      rownames(lmr)<- NULL
   ier<- apply(mdlB, 2, mean, na.rm=TRUE)
      ier<- data.frame(id=names(ier), r=ier)
      rownames(ier)<- NULL
   ler<- apply(mdlB, 1, mean, na.rm=TRUE)
      st<- strsplit(names(ler), "_")
      st<- data.frame(chr=sapply(st, function(x)x[1]), pos=sapply(st,function(x)as.integer(x[2])))
      ler<- cbind(st, r=ler)
      rownames(ler)<- NULL

   fstr<- paste(bcftools, " +fill-tags -r ", chr, sep="")
      fstr<- paste(fstr, " ", fl, sep="")
      fstr<- paste(fstr, " -- -t MAF ", sep="")
      fstr<- paste(fstr, " | ", bcftools, " query -f '%POS\t%INFO/MAF\n'", sep="")
   maf<- system(fstr, intern=TRUE)
      maf<- strsplit(maf, "\t")
      maf<- do.call("rbind", maf)
      maf<- as.data.frame(maf)
      colnames(maf)<- c("pos", "MAF")
      maf$pos<- as.integer(maf$pos)
      maf$MAF<- as.numeric(maf$MAF)
   idx<- match(lmr$pos, maf$pos)
   lmr$f<- maf$MAF[idx]
   
   list(imr=imr, lmr=lmr, ier=ier, ler=ler)

   # Create output directory if it doesn't exist
   if (!dir.exists(outdir)) {
      dir.create(outdir)
   }

   # Export each element to a CSV file
   write.csv(imr, file.path(outdir, paste("imr", chr, sep = "_")), row.names = FALSE)
   write.csv(lmr, file.path(outdir, paste("lmr", chr, sep = "_")), row.names = FALSE)
   write.csv(ier, file.path(outdir, paste("ier", chr, sep = "_")), row.names = FALSE)
   write.csv(ler, file.path(outdir, paste("ler", chr, sep = "_")), row.names = FALSE)

}

#######################
## Usage
#######################

# Read in arguments
args <- commandArgs(TRUE) 
fl <- args[1]
pedigree <- args[2]
chr <- args[3]
tmpdir <- args[4]
outdir <- args[5]
bcftools <- args[6]

ped<- read.csv(pedigree, header=TRUE, stringsAsFactors=FALSE)
colnames(ped)[1]<- "id"


mdlf(fl=fl, ped=ped, outdir=outdir, chr=chr, tmpdir=tmpdir)
