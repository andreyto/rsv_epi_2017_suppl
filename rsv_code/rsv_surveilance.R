## Utils

replace_file_ext <- function(x,new_ext) {
  file.path(dirname(x),paste0(sub("\\.[^.]*$","",basename(x)),new_ext))
}

rsv.vac.debug.mode = F

id.from.seq.names <- function(defline) {
  library(stringr)
  str_extract(defline,"^[^\\s]+")
}

sanitize.seq.id.for.phangorn <- function(id) {
  gsub("[-|]","_",id)
}

rsv.genotype.from.id <- function(id) {
  library(tidyr)
  tidyr::extract(data.frame(names=names(s)),"names","id",'^[^-]+-(\\S*) ')
}

seq.dist.to.ident <- function(d,dist.transform=c("sqrt","none")) {
  dist.transform = dist.transform[[1]]
  if(inherits(d,"dist")) {
    d = as.matrix(d)
  }
  if(dist.transform=="sqrt") {
    1 - d**2
  }
  else {
    1 - d
  }
}

ident.to.seq.dist <- function(i,dist.transform=c("sqrt","none")) {
  dist.transform = dist.transform[[1]]
  if(dist.transform=="sqrt") {
    sqrt(1 - i)
  }
  else {
    1 - i
  }
}

## Pick exemplar for each group produced by cluster.with.dist.cutoff()
## Criteria is min number of degenerates, then min distance to other members in the cluster
pick.dist.group.exemplar <- function(d,group,seq.attr) {
  library(data.table)
  dgroup = data.table(rn=names(group),group=group)
  setkey(dgroup,rn)
  if(is.null(seq.attr$Count)) {
    seq.attr[,Count:=1]
  }
  if(is.null(seq.attr$n_nongenic)) {
    seq.attr[,n_nongenic:=0]
  }
  
  dind = seq.attr[,.(rn,n_nongenic,Count)]
  setkey(dind,rn)
  
  if(!is.null(d)) {
    dind = dind[labels(d),nomatch=0]
    stopifnot(nrow(dind)==length(labels(d)))
    dm = as.matrix(d)
  }
  else {
    dm = NULL
  }
  
  dind[,ind:=.I]
  dgroup = dgroup[dind,on="rn"]
  setkey(dgroup,group)
  ## below, we compute, for each sequence in a group (row),
  ## the total sum of distances to each other sequence in the same group.
  ## That would be done by extracting a block from the distance matrix dm[ind,ind],
  ## and taking rowSums.
  ## We also need to multiply each distance by the size of the group that was previously
  ## assigned to each sequence (weigh by size), by muliplying each column of a distance block
  ## by the n_group of the column. That would be done as t(t(m)*v). We remove the outer t() by
  ## replacing the final rowSums with colSums, and, becaus the m is symmetrical, there is no
  ## need to do the inner t().
  ## If the distance matrix is not provided, we assume zero distance for all pairwise comparisons
  dist_group = dgroup[,.(ind,
                         distSum=if(!is.null(dm)) colSums(dm[ind,ind,drop=F]*dgroup[ind,Count]) else 0),
                      by=group]
  dgroup = dgroup[dist_group[,.(ind,distSum)],on="ind"]
  dex = dgroup[,
               .(ind=ind[order(distSum,n_nongenic)][1],
                 members=.N,
                 members_with_nongenic=sum(n_nongenic>0),
                 Count=sum(Count)),
               by=group]
  
  dex$rn_ex = dind[dex$ind,rn]
  dind = dind[dgroup,on="rn"][,.(rn,group)]
  dind = dind[dex,on="group"]
  
  dind_ex = dind[dex$ind]
  dind_ex$members = dex$members
  dind_ex$members_with_nongenic = dex$members_with_nongenic
  dind_ex$Count =dex$Count
  
  list(dind=dind,dind_ex=dind_ex)
}


duplicated_index <- function(x) {
  library(data.table)
  xord = order(x)
  x = x[xord]
  ord = data.table::data.table(iid=seq(length(xord)))[xord]
  y = data.table::data.table(is_uni=!duplicated(x),ord)
  y[,idgr_uniq:=cumsum(is_uni)]
  y_uniq = y[is_uni==T,.(idgr_uniq,iid_uniq=iid)]
  y_uniq[y,on="idgr_uniq"]
}

duplicated_index_names <- function(x) {
  dup_ind = duplicated_index(x)
  nm = names(x)[dup_ind[,iid]]
  group = dup_ind$idgr_uniq
  names(group) = nm
  group
}

init_seq_for_dedup <- function(seq) {
  attr = seq$attr
  attr[,rn:=SeqID]
  attr[,Count:=1]
  attr_all = copy(attr)
  attr_all[,rn_ex:=rn]
  seq$attr = attr
  seq$attr_all = attr_all
  seq
}


dedup_seq <- function(seq,attr_aggr_fun) {
  attr = seq$attr
  s = seq$seq
  stopifnot(!anyDuplicated(names(s)))
  group = duplicated_index_names(s)
  ex = pick.dist.group.exemplar(d=NULL,group=group,seq.attr=attr)
  aggregate_seq(seq=seq,exemplars = ex$dind,attr_aggr_fun = attr_aggr_fun)
}

aggregate_seq <- function(seq,exemplars,attr_aggr_fun) {
  stopifnot(all(sort(unique(seq$attr_all$rn_ex))==sort(exemplars$rn)))
  attr = copy(seq$attr)
  s = seq$seq
  stopifnot(!anyDuplicated(names(s)))
  attr_all = copy(seq$attr_all)
  ## for sanity check at the end of this function
  total = attr_all[,sum(Count)]
  ##GOTCHA in data.table. Use the two lines below to rename a field.
  ##If instead you use setnames(attr_all,"rn_ex","rn_ex_old"), and
  ##"rn_ex_old" already exists, data.table will simply delete the "rn_ex"
  ##and use the current values in "rn_ex_old" (2017-11-21)
  attr_all[,rn_ex_old:=rn_ex]
  attr_all[,rn_ex:=NULL]
  attr_all = attr_all[exemplars[,.(rn_ex,rn)],on=c(rn_ex_old="rn")]
  attr_all[,SeqID_uniq:=rn_ex]
  gr_summ = attr_all[,.(Count=sum(Count),SeqID_all=paste(SeqID,collapse = ",")),by=rn_ex]
  attr = gr_summ[attr,on=c(rn_ex="rn"),nomatch=0]
  attr[,rn:=rn_ex]
  attr[,rn_ex:=NULL]
  s = s[attr$rn]
  seq = do.call(attr_aggr_fun,
                list(list(seq=s,attr=attr,attr_all=attr_all))
  )
  stopifnot(all(sort(seq$attr$rn)==sort(unique(exemplars$rn_ex))))
  stopifnot(seq$attr_all[,sum(Count)]==total)
  stopifnot(seq$attr[,sum(Count)]==total)
  seq
}

rsv_aggr_attr_exemplar <- function(seq) {
  attr = seq$attr
  s = seq$seq
  attr_all = seq$attr_all
  gr_summ = rsv_aggr_attr(attr_all,by="rn_ex")
  attr = gr_summ[attr,on=c(rn_ex="rn"),nomatch=0]
  attr[,rn:=rn_ex]
  attr[,rn_ex:=NULL]
  s = s[attr$rn]
  stopifnot(all(sort(attr$rn)==sort(gr_summ$rn_ex)))
  list(seq=s,attr=attr,attr_all=attr_all)
}

rsv_aggr_attr <- function(attr,by) {
  gr_summ = attr[,.(Hospital_Stay_Long_Ratio=mean(Hospital_Stay_Long==">24H",na.rm=T),
                    Age_month_median=median(Age_month,na.rm=T),
                    Count=sum(Count)),by=by]
  gr_summ[,Age_group_median:=cut(Age_month_median,
                                 c(0,12,65*12,1e6),
                                 include.lowest = F,
                                 ordered_result = T,
                                 labels=c("Infant","ChildOrAdult","Senior"))]
  gr_summ[,Age_month_qua:=quantcut.ordered(Age_month_median)]
  gr_summ
}

exemplar.report <- function(exemplars) {
  exemplar = exemplars$dind_ex
  report.section = report$add.header("Effect of the clustering",section.action="push",sub=T)
  y = exemplar[,.(has_nongenic=members_with_nongenic>0,singleton=members==1)]
  t = xtabs(~has_nongenic+singleton,exemplar[,.(has_nongenic=members_with_nongenic>0,singleton=members==1)])
  report$add.printed(t,caption="Count of sequences with nongenic symbols in singletons vs larger clusters")
  library(DescTools)
  report$add.printed(DescTools::PercTable(t,margins = 1:2,rfrq = "111"),
                     caption = "Proportion of sequences with nongenic symbols in singletons vs larger clusters")
  report$add.printed(summary(t),caption = "Significance of the difference")
  report$pop.section()
}

alignment.as.bios2mds <- function (seq, aa.to.upper = TRUE) 
{
  if(!inherits(seq,"character")) seq = as.character(seq)
  rn = names(seq)
  if (aa.to.upper) 
    seq <- toupper(seq)
  seq <- strsplit(seq, split = "")
  names(seq) = rn
  class(seq) <- c("align")
  return(seq)
}

## Compute distance matrix from alignment.
## This respects all masks set on alignment.
## Note: bios2mds.sim only works with 20 standard codes plus gap
## Note: seqinr.sim uses Fitch matrix (the number of nucelotide changes required to change one AA into another)
## something might be odd with seqinr::dist.alignment(ali.ape,matrix = "similarity") -
## it returns zero distances for sequences with differences. TODO: check that this is not due to ignoring
## pairwise positions with degenerate symbols. With DECIPHER, this is indeed the case. That makes
## the dissimilarity to to be not a distanceby failing the triangle enequality - it is possible to
## have AB=0 while AC!=BC.

rsv.ali.deriv <-function(ali,
                         type=c("DECIPHER.sim","seqinr.ident","seqinr.sim","bios2mds.sim","phangorn.sim","DECIPHER.ident"),
                         phangorn.model="JTT",
                         bios2mds.sub.mat.id="PAM250",
                         mask.gaps=T,
                         mask.gaps.min.fraction=0.5,
                         mask.gaps.min.block.width=4) {
  dist.transform = "none"
  ##Note that removing gaps is needed to build
  ##good trees as per msa package advice, but might show zero distance
  ##for different sequences in plots and tables
  type = type[[1]]
  if(mask.gaps) {
    ali = mask.gaps.ali(ali, 
                        min.fraction=mask.gaps.min.fraction, 
                        min.block.width=mask.gaps.min.block.width)
  }
  if(type %in% c("seqinr.ident","seqinr.sim")) {
    dist.transform = "sqrt"
    library(seqinr)
    matrix = switch(type,
                    seqinr.ident="identity",
                    seqinr.sim="similarity",
                    stop("Unknown similarity type"))
    #ali.m = as.matrix(ali.masked)
    #ali.ape = ape::as.alignment(ali.m)
    ## this removes everything that was masked
    ali.s = as.character(ali)
    ali.ape = seqinr::as.alignment(nb=length(ali.s),
                                   nam=names(ali.s),
                                   seq=ali.s)
    ##With Bioc 3.3 we can instead use:
    #msa::msaConvert(ali, type="seqinr::alignment")
    d = seqinr::dist.alignment(ali.ape,matrix = matrix)
  }
  else if(type=="phangorn.sim") {
    #ali = masked.copy.ali(ali,drop.rows = T,drop.columns = T,deep = T)
    phdat = alignment.as.phyDat(ali)
    d = phangorn::dist.ml(phdat,model=phangorn.model)
  }
  else if(type=="bios2mds.sim") {
    library(bios2mds)
    ## this removes everything that was masked
    ali.s = as.character(ali)
    d = as.dist(bios2mds::mat.dis(ali.s,ali.s,sub.mat.id = bios2mds.sub.mat.id,sqrt=T))
    stopifnot(all(names(ali.s)==labels(d)))
  }
  else if(type %in% c("DECIPHER.sim","DECIPHER.ident")) {
    library(DECIPHER)
    correction = switch(type,
                        DECIPHER.ident="none",
                        DECIPHER.sim="Jukes-Cantor", #for AA, this is actually called Bishop–Friday model,
                        #method's code is aware of AAStringSet. Correction for multiple substiutions makes
                        #distances larger, and the corrected output might have value > 1.0
                        stop("Unknown similarity type"))
    ## this removes everything that was masked
    ali.s = as(ali,"XStringSet")
    d = as.dist(DECIPHER::DistanceMatrix(ali.s,correction = correction))
    stopifnot(all(names(ali.s)==labels(d)))
  }
  else {
    stop(sprintf("Unknown value of type: %s",type))
  }
  list(d=d,dist.transform=dist.transform)
}

#' this will always construct cophenetic distance and return it as dist component of returned list object
rsv.make.tree.and.dist <- function(dist,
                                   method=c("nj","upgma","load"),
                                   root.tree=T,
                                   tree.file=NULL) {
  method = method[[1]]
  #t = ape::bionj(dist)
  if(method=="nj") {
    #library(ape)
    #t = ape::nj(dist)
    library(phangorn)
    t = phangorn::NJ(dist)
  }
  else if(method=="upgma") {
    library(phangorn)
    t = phangorn::upgma(dist) #this just calls hclust
  }
  else if(method=="load") {
    if(is.null(tree.file)) {
      stop("Method is 'load' but tree.file is NULL")
    }
    t = load.tree.raxml(tree.file)
  }
  else stop(sprintf("Unknown value of method: %s"),method)
  #t = fastme.bal(d)
  #SeqGroup = seq$attr[labels(ali.der$d),"SeqGroup"]
  #outgroup = labels(ali.der$d)[SeqGroup!=target.SeqGroup]
  #group = labels(dist)[grepl("_[BC]$",labels(dist))]
  #node.med = ape::getMRCA(t,group)
  node.med = get.tree.node.medoid(t,format="index")
  if(root.tree) try({
    #t = ape::root(t,node=node.med, r=T)
    t = ape::root(t,outgroup=node.med, r=T)
  })
  list(dist = as.dist(ape::cophenetic.phylo(t)),
               dist.transform = "sqrt",
               tree = t,
               node.med = node.med)
}

rsv.ali.to.tree.attr <- function(ali,
                                 seq.attr,
                                 target.SeqGroup=NULL,
                                 make.tree=T,
                                 use.tree.dist=F,
                                 tree.method="nj",
                                 root.tree=T,
                                 mask.gaps=T,
                                 mask.gaps.min.fraction=0.5,
                                 mask.gaps.min.block.width=4,
                                 tree.file=NULL,
                                 tree.dist.type="DECIPHER.sim",
                                 ident.dist.type="DECIPHER.ident") {
  seq.attr = copy(seq.attr)
  if(mask.gaps) {
    ali = mask.gaps.ali(ali, 
                        min.fraction=mask.gaps.min.fraction, 
                        min.block.width=mask.gaps.min.block.width)
  }
  ##gaps are already masked above
  ali.der = rsv.ali.deriv(ali,type=tree.dist.type,mask.gaps = F)
  dist.transform = ali.der$dist.transform
  if(use.tree.dist) make.tree = T
  if(make.tree) {
    tree.dist = rsv.make.tree.and.dist(ali.der$d,
                                       method=tree.method,
                                       root.tree = root.tree,
                                       tree.file = tree.file)
  }
  else {
    tree.dist = NULL
  }
  ali.der.ident = rsv.ali.deriv(ali,type=ident.dist.type,mask.gaps = F)
  ident.mat = as.matrix(
    ##gaps are already masked above
    seq.dist.to.ident(ali.der.ident$d,
                      dist.transform = ali.der.ident$dist.transform)
  )[seq.attr$rn,seq.attr$rn]
  
  ##can instead use cophenetic distance tree.dist$dist (note that it can somehow be negative
  ##when consensus sequence is added to the alignment)
  d.mat = as.matrix(
    if(use.tree.dist) tree.dist$dist else ali.der$d
  )[seq.attr$rn,seq.attr$rn]
  if(use.tree.dist) dist.transform = tree.dist$dist.transform
  ## mean dist from each sequence to all other sequences, for plotting
  ## d.mat is zero at the diagonal, so we can sum the entire row and divide 
  ## by the number of "other" elements. But we subtract the diagonal
  ## for both distance and identity cases for generality, although it
  ## is only needed for the identity.
  seq.attr$d.mean = (rowSums( d.mat ) - diag(d.mat))/(ncol(d.mat)-1)
  seq.attr$ident.mean = (rowSums( ident.mat ) - diag(ident.mat))/(ncol(ident.mat)-1)
  
  #TODO: fix for self-identity and propagate this and code above to HBV
  if(!is.null(target.SeqGroup)) {  
    ## mean dist from each sequence to sequences in the target group, for plotting
    d.mat.all.to.sel = d.mat[,seq.attr$SeqGroup==target.SeqGroup]
    seq.attr$d.mean.sel = rowSums( d.mat.all.to.sel )/(ncol(d.mat.all.to.sel)-1)
    
    ident.mat.all.to.sel = ident.mat[,seq.attr$SeqGroup==target.SeqGroup]
    seq.attr$ident.mean.sel = rowSums( ident.mat.all.to.sel )/(ncol(ident.mat.all.to.sel)-1)
  }
  
  return (list(tree.dist=tree.dist,seq.attr=seq.attr,
               d.mat=d.mat,ident.mat=ident.mat,dist.transform=dist.transform))
}

rsv.select.medoid <- function(d.mat,seq.attr,target.SeqGroup) {
  stopifnot(all(rownames(d.mat)==seq.attr$rn))
  filt.col = (seq.attr$SeqGroup==target.SeqGroup)
  filt.row = (filt.col & (seq.attr$n_nongenic==0))
  d.mat.sel = d.mat[filt.row,filt.col]
  d.mean.sel = rowSums( d.mat.sel )/(ncol(d.mat.sel)-1)
  #d.mean.sel = robustbase::rowMedians(d.mat.sel) ##keeps names
  ##medoid
  name.med = names(which.min(d.mean.sel))[1]
  return (list(name.sel=name.med))
}

rsv.select.nearest <- function(d.mat,seq.attr,target.SeqGroup,name.sel) {
  stopifnot(all(rownames(d.mat)==seq.attr$rn))
  filt.row = (
    (seq.attr$SeqGroup==target.SeqGroup) & 
      (seq.attr$n_nongenic==0) &
      seq.attr$rn != name.sel
  )
  
  d.mat.sel = d.mat[filt.row,name.sel]
  name.sel.targ = name.sel
  min.ind = which.min(d.mat.sel)[1]
  name.sel = names(min.ind)
  return (list(name.sel=name.sel, name.sel.targ=name.sel.targ, d.sel.targ=d.mat.sel[min.ind]))
}


rsv.make.consensus <- function(ali,seq.attr,target.SeqGroup,update.ali=T,seq.attr.copy.fields=c("Genotype","SeqGroup")) {
  setkey(seq.attr,"rn")
  seq.attr = seq.attr[rownames(ali)]
  stopifnot(all(rownames(ali)==seq.attr$rn))
  row.sel=(seq.attr$SeqGroup==target.SeqGroup)
  if(length(row.sel)==0) {
    return (list(seq.sel=NULL))
  }
  else {
    rmask = rowmask(ali)
    cmask = colmask(ali)
    rowmask(ali,append="union") = IRanges(start=!row.sel)
    ##we do not want to get any "no consensus" letters(e.g. `?`) in the output, so 
    ##trying to set threshold but it has strange effects - both low and high values
    ##result in lots of `?`. Using our wrapper that takes majority vote at threshold=0.
    ##TODO: read the code of consensusString to see how
    ##threshold is used.
    seq.sel = get.consensus.from.ali(ali,drop.rows=T,drop.columns=F,threshold=0)
    name.sel = sprintf("Consensus_%s",target.SeqGroup)
    names(seq.sel) = name.sel
    seq.sel = call.ctor(unmasked(ali),
                        seq.sel,
                        use.names=T)
    if(!update.ali) {
      return (list(seq.sel=seq.sel))
    }
    else {
      seqs = unmasked(ali)
      ali = call.ctor(ali,
                      append(seqs,
                             seq.sel),
                      rowmask=rmask,
                      colmask=cmask
      )
      attr.sel = copy(seq.attr[which(row.sel)[1]])
      ind.zap.fields = which(!(colnames(attr.sel) %in% seq.attr.copy.fields))
      for(ind in ind.zap.fields) attr.sel[,ind] = NA
      attr.sel[,rn:=name.sel]
      seq.attr = rbind(seq.attr,attr.sel)
      setkey(seq.attr,"rn")
      return (list(name.sel=name.sel,ali=ali,seq.attr=seq.attr))
    }
  }
}

rsv.make.internal.node <- function(ali,seq.attr,
                                   target.SeqGroup,
                                   tree.target.SeqGroup.only=F,
                                   update.ali=T,
                                   seq.attr.copy.fields=c("Genotype","SeqGroup")) {
  setkey(seq.attr,"rn")
  seq.attr = seq.attr[rownames(ali)]
  stopifnot(all(rownames(ali)==seq.attr$rn))
  row.sel=(seq.attr$SeqGroup==target.SeqGroup)
  if(length(row.sel)==0) {
    return (list(seq.sel=NULL))
  }
  else {
    rmask = rowmask(ali)
    cmask = colmask(ali)
    if(tree.target.SeqGroup.only) {
      rowmask(ali,append="union") = IRanges(start=!row.sel)
    }
    ali.targ = masked.copy.ali(ali,drop.rows = T,drop.columns = F)
    seq.attr.targ = filter.seq.attr.by.ali(ali.targ,seq.attr,drop.rows=T)
    tree.attr = rsv.ali.to.tree.attr(ali=ali.targ,
                                     seq.attr = seq.attr.targ,
                                     target.SeqGroup=(
                                       if(tree.target.SeqGroup.only) 
                                         NULL else target.SeqGroup),
                                     make.tree=T,
                                     use.tree.dist=F,
                                     tree.method="nj",
                                     root.tree=F)
    anc = reconstruct.ancestors(ali.targ,tree.attr$tree.dist$tree)
    #node.sel = get.tree.mrca(anc$tree,
    #              tips=seq.attr.targ[SeqGroup==target.SeqGroup]$rn,
    #              format="label")
    
    ##This takes center of entire tree
    cot.type = "full_tree"
    if (cot.type == "full_tree") {
      node.sel = get.tree.node.medoid(anc$tree,
                                      format="label")
    }
    ##This takes center of the target group
    else {
      node.sel = get.center.of.tree(anc$tree,
                                    nodes=seq.attr.targ[SeqGroup==target.SeqGroup]$rn,
                                    format="label")
    }
    seq.sel = anc$seq[node.sel]
    name.sel = sprintf("CenterOfTree_%s",target.SeqGroup)
    names(seq.sel) = name.sel
    seq.sel = call.ctor(unmasked(ali),
                        seq.sel,
                        use.names=T)
    if(!update.ali) {
      return (list(seq.sel=seq.sel))
    }
    else {
      seqs = unmasked(ali)
      ali = call.ctor(ali,
                      append(seqs,
                             seq.sel),
                      rowmask=rmask,
                      colmask=cmask
      )
      attr.sel = copy(seq.attr[which(row.sel)[1]])
      ind.zap.fields = which(!(colnames(attr.sel) %in% seq.attr.copy.fields))
      for(ind in ind.zap.fields) attr.sel[,ind] = NA
      attr.sel[,rn:=name.sel]
      seq.attr = rbind(seq.attr,attr.sel)
      setkey(seq.attr,"rn")
      return (list(name.sel=name.sel,ali=ali,seq.attr=seq.attr))
    }
  }
}

rsv.selection.stats <- function(d.mat,ident.mat,name.sel) {
  stopifnot(all(rownames(d.mat)==rownames(ident.mat)))
  df = data.table(rn=rownames(d.mat))
  df$d.sel = d.mat[,name.sel]
  df$ident.sel = ident.mat[,name.sel]
  ## exclude self
  ##TODO: propagate to HBV
  df = df[!(rn %in% name.sel)]
  setkey(df,"rn")
  df
}

seq_letter_frequency <- function(seq,standard_alphabet=Biostrings::AA_STANDARD){
  lt.freq = colMeans(letterFrequency(seq, as.prob=T, letters = setdiff(alphabet(seq),"-")))
  lt.freq = data.table(
    Letter=factor(names(lt.freq),levels=names(lt.freq)),
    Frequency=lt.freq)
  nz_lett = as.character(lt.freq[Frequency>0,Letter])
  not_in_seq = setdiff(standard_alphabet,nz_lett)
  not_in_standard_alphabet = setdiff(nz_lett,standard_alphabet)
  list(freq=lt.freq,
       nz_lett=nz_lett,
       seq_and_standard_alphabet=intersect(nz_lett,standard_alphabet),
       not_in_seq=not_in_seq,
       not_in_standard_alphabet=not_in_standard_alphabet)
}

rsv.sequence.qc.report <- function(seq,seq.attr,descr=NULL,standard_alphabet = Biostrings::AA_STANDARD) {
  stopifnot(length(seq)==nrow(seq.attr))
  caption = "QC of sequences"
  if(!is.null(descr)){
    caption = paste(caption,descr)
  }
  report.section = report$add.header(caption,section.action="push",sub=T)
  fr = seq_letter_frequency(seq)
  lt.freq = fr$freq
  report$add.table(lt.freq,caption="Letter frequency")
  report$add.descr(sprintf("Used letters: %s",paste(fr$nz_lett,collapse = ",")))
  report$add.descr(sprintf("%s letters out of %s in the standard AA alphabet are used, 
                           the following letters are not used: %s,
                           the following letters are used but are not in the standard alphabet: %s",
                           length(fr$seq_and_standard_alphabet),
                           length(standard_alphabet),
                           paste(fr$not_in_seq,collapse = ","),
                           paste(fr$not_in_standard_alphabet,collapse = ",")))
  report$add(ggplot(aes(x=Letter,y=Frequency),
                    data = lt.freq) + 
               geom_bar(stat="identity"),
             caption = "Letter frequency in sequences"
  )
  tot_nongenic = sum(seq.attr$n_nongenic>0)
  report$add.descr(sprintf("**Total count of sequences with degenerate and other special characters present is %s,
                           which is %g%% of all %s sequences**",tot_nongenic,(tot_nongenic/nrow(seq.attr)*100),nrow(seq.attr)))
  
  report$pop.section()
}




load.tree.raxml <- function(file.name) {
  ggtree::get.tree(ggtree::read.raxml(file.name))
}

plot.rsv.tree.raxml <- function() {
  caption = sprintf("RSV tree properties")
  report.section = report$add.header(
    caption,
    section.action = "push")
  
  library(ggtree)
  
  tree = read.raxml("RAxML_bipartitionsBranchLabels.TGAMMAHIVWF")
  state = str_sub(str_split_fixed(tree@phylo$tip.label,"_",3)[,2],1,2)
  tree = groupOTU(tree,split(tree@phylo$tip.label,state))
  pl = ggtree(tree,aes(color=group),layout='circular') + geom_label(aes(label=bootstrap, fill=bootstrap)) + geom_tiplab() +
    scale_fill_continuous(low='darkgreen', high='red') + theme_tree2(legend.position='right')
  
  print(pl)
}

rsv_nuc_to_prot <- function(nuc_files_patt,nuc_ali_dir,
                            msa.method = "MAFFT",
                            msa.args = list(param="--maxiterate 1000 --localpair --thread -1")) {
  inp_fas = list.files(dirname(nuc_files_patt),basename(nuc_files_patt),full.names = T)
  dir.create(nuc_ali_dir)
  for(inp_fa in inp_fas) {
    out_ali = replace_file_ext(file.path(nuc_ali_dir,basename(inp_fa)),".ali.fasta")
    seq = Biostrings::readDNAStringSet(inp_fa)
    ali = do.call(mgsat.msa,
                  c(list(seq,msa.method=msa.method),msa.args))
    Biostrings::writeXStringSet(as(ali,"XStringSet"),out_ali,format = "fasta")
  }
}

rsv_pdb_radar_ref_ali <- function(pdb_seq,ref_seq,pdb_out_dir) {
  pdb_base = replace_file_ext(basename(pdb_seq),"")
  out_dir = file.path(pdb_out_dir,pdb_base)
  dir.create(out_dir,recursive = T)
  out_ali = replace_file_ext(file.path(out_dir,paste0("pdb_",basename(ref_seq))),".ali.fasta")
  msa.method = "MAFFT"
  msa.args = list(param="--maxiterate 1000 --localpair")  
  s_pdb = Biostrings::readAAStringSet(pdb_seq)
  seq = c(s_pdb,
          Biostrings::readAAStringSet(ref_seq)
  )
  ali = do.call(mgsat.msa,
                c(list(seq,msa.method=msa.method),msa.args))
  Biostrings::writeXStringSet(as(ali,"XStringSet"),out_ali,format = "fasta")
  return (out_ali)
}

ordered.color.legend.report <- function(tbl_col) {
  title = "Variation Frequency"
  pal = ordered.color.legend(tbl_col,title=title)
  caption=sprintf("Palette for %s",title)
  report$add.table(pal$tbl_palette,caption = caption)
  report$add(pal$plot_palette,caption = caption)
}

read_RADAR_var_report <- function(file_name) {
  library(data.table)
  vars = data.table::fread(file_name,colClasses = "character")
  #skip footer
  vars = vars[1:(which(str_blank(vars$Name))[1]-1)]
  ## sanity check that first row looks as expected
  stopifnot(vars[1,grepl("^Ref",Name)])
  var_cols = grepl("^[A-Z]+_[0-9]+",colnames(vars))
  ref_freq = vars[Name=="Ref Frequency %"]
  stopifnot(nrow(ref_freq)<=1)
  if(nrow(ref_freq)>0) {
    stopifnot(vars[1:3,grepl("^Ref",Name)])
    
    ref_freq = ref_freq[,var_cols,with=F]
    stopifnot(all(grepl("%$",ref_freq)))
    ref_freq_names = colnames(ref_freq)
    ref_freq = as.numeric(sub("%$","",ref_freq))
    names(ref_freq) = ref_freq_names
    ##TODO: compute from absolute counts in the previous string in case RADAR will round up very small % to 100
    ref_freq = 1 - (ref_freq / 100)
    vars_bases_rows = -(1:3)
  }
  else {
    stopifnot(!any(vars[2:3,grepl("^Ref",Name)]))
    ref_freq = NULL
    vars_bases_rows = -1
  }
  vars_bases = vars[vars_bases_rows,]
  bases = vars_bases[,var_cols,with=F]
  colnames_bases = colnames(bases)
  rownames_bases = rownames(bases)
  bases = as.matrix(bases)
  bases = matrix(stringr::str_trim(bases),nrow=nrow(bases),ncol=ncol(bases))
  rownames(bases) = vars_bases[,Name]
  colnames(bases) = colnames_bases
  list(bases=bases,vars_freq=ref_freq)
}

#  vars_tbl = data.table(pos=vars_pos,
#freq=vars_freq,
#name=names(vars_freq))

plot.variation.freq <- function(vars_tbl,pos=c("linear","compressed")) {
  pos = pos[1]
  pl.data = copy(vars_tbl)
  if(pos=="compressed") {
    pl.data = pl.data[freq>0][,pos:=ordered(pos)]
  }
  pl.data[,freq:=freq*100]
  min.x = min(pl.data$freq)
  max.x = max(pl.data$freq)
  pl.data[,min.x:=min.x]
  pl.data[,max.x:=max.x]
  pl = ggplot(pl.data, aes(x = freq, y = pos)) +  
    geom_point(size=rel(4)) +
    geom_segment(aes(xend = freq,yend=pos,y=pos,
                     x=0)) +
    #geom_bar(stat="identity",width=0.8) +
    scale_x_continuous(limits=c(
      min.x-(max.x-min.x)*0.1,
      max.x)
    ) +
    xlab("Variation frequency (%)") +
    ylab("Residue position") +
    #coord_flip() +
    theme(plot.title = element_text(size = rel(2)),
          axis.text.x = element_text(color=c("black","black"),size = rel(1)), #,angle=90
          axis.text.y = element_text(size = rel(1)))
  
  return(pl)
}

rsv_pdb_selections_RADAR <- function(radar_var,pdb_out_dir,freq_type=c("RADAR","REF","CONS")) {
  ## Extract variable, non-gap, non-degen positions from RADAR variant file,
  ## and save them. These positions are used without remapping on a PDB
  ## file that had its residue numbering remapped with the output from
  ## rsv_ali_to_index(). Such approach is only valid when both A and B
  ## references have no indels relative to each other.
  radar = read_RADAR_var_report(radar_var)
  bases = radar$bases
  ## mark all degenerates as NA (more than one letter per base or X) and "sequencing"
  ## (as opposed to alignment) gaps (represented in RADAR as dots). The dot treatment
  ## might be specific to RSV F where such things mean sequencing
  ## errors or truncated sequences, but is likely to be generic as well.
  count_alleles = matrix(stringr::str_length(bases),nrow=nrow(bases),ncol=ncol(bases))
  bases[which(count_alleles>1,arr.ind = T)] = NA
  bases[which(bases==".",arr.ind = T)] = NA
  bases[which(bases=="X",arr.ind = T)] = NA
  if( freq_type == "FREQ") {
    vars_freq = colMeans(bases!="",na.rm = T)
    ##we should not have any columns where everything is NA (why not?)
    stopifnot(all(!is.na(vars_freq)))
  }
  else if(freq_type == "CONS") {
    #x = melt(as.data.table(bases,keep.rownames = T),id.vars = "rn",variable.name = "pos",value.name = "base")
    vars_freq = apply(bases,2,function(x) { 
      x = data.table(base=x) 
      x = x[!is.na(base),.(cnt=.N),by=base]
      1- x[,max(cnt)]/x[,sum(cnt)] 
    })
  }
  else if( freq_type == "RADAR") {
    vars_freq = radar$vars_freq
  }
  else {
    stop(sprintf("Unknown freq_type argument: %s",freq_type))
  }
  vars_pos = names(vars_freq)
  vars_pos = as.numeric(stringr::str_match(vars_pos,"_([0-9]+)$")[,2])
  vars_tbl = data.table(pos=vars_pos,
                        freq=vars_freq,
                        name=names(vars_freq))
  if(vars_tbl[freq==0,.N]==0) {
    stop("RADAR variant table contains no 100% conserved positions - suspecting that you are
         loading a compressed (filtered to variations) table. We need the unfiltered table here.")
  }
  is_cons = (bases=="")
  cnt_var = rowSums(!is_cons,na.rm = T)
  ref_len = rowSums(!is.na(is_cons))
  var_rows = data.table(cnt_var_rank=rank(cnt_var,ties.method="first"),name=names(cnt_var),cnt_var=cnt_var,ref_len=ref_len)
  var_rows = var_rows[order(cnt_var_rank,decreasing = T)]
  var_rows[,n_seq:=.N]
  var_rows[,cnt_var_rank_pct:=cnt_var_rank/n_seq*100]
  var_rows[,cons_pct:=(ref_len-cnt_var)/ref_len*100]
  cons_tot_pct = var_rows[,mean(cons_pct)]
  cum_var = var_rows[,.(rank_pct=max(cnt_var_rank_pct)),by=cnt_var]
  ##shift cnt_var down; then rank_pct means that many records have less than cnt_var
  cum_var[,cnt_var:=c(cnt_var[1]+1,cnt_var[1:(.N-1)])]
  report$add.table(var_rows,caption = sprintf("Ranked per-sequence conservation. Mean conservation is %s%%",cons_tot_pct))
  report$add.table(cum_var,caption = "Cumulative percentage ranking of the number of variations per sequence. 
                   Percentage value means that many sequences have number of variations strictly less than cnt_var.")
  out_vars_pos = replace_file_ext(file.path(pdb_out_dir,paste0("vars_pos_",basename(radar_var))),".txt")
  report$add(plot.variation.freq(vars_tbl,pos="linear"),
             caption = "Variation frequency linear coordinates",
             width=800,
             height=1000)
  report$add(plot.variation.freq(vars_tbl,pos="compressed"),
             caption = "Variation frequency compressed coordinates",
             width=800,
             height=1000)
  list(out_vars_pos=out_vars_pos,vars_tbl=vars_tbl)
}


rsv_pdb_selections <- function(input_type=c("RADAR","MSA"),
                               radar_var=NULL,
                               msa=NULL,
                               pdb_out_dir,
                               freq_type="RADAR") {
  input_type = input_type[1]
  if(input_type=="RADAR") {
    res_var = rsv_pdb_selections_RADAR(radar_var=radar_var,pdb_out_dir=pdb_out_dir,freq_type=freq_type)
  }
  else {
    stop(sprintf("Not implemented for input_type: %s",input_type))
  }
  write.table(res_var$vars_tbl[freq>0],res_var$out_vars_pos,row.names = F,sep="\t")
  res_var$out_vars_pos
}

#' adhoc analysis for epi manuscript 2017-12-01 report
rsv_adhoc_get_top_cumulative_clusters_weight <- function(x,n_top=5) {
  y = x[,.(Count=sum(Count)),by=SeqID][order(Count,decreasing = T)]
  y[,Subtype:=substring(SeqID,stringr::str_length(SeqID))]
  get_top_ratio = function(z) {z=sort(z,decreasing = T); sum(z[1:n_top])/sum(z)}
  z = y[,.(TopRatio=get_top_ratio(Count)),by=Subtype]
  z
}

rsv_make_vars_colors <- function(vars_files,out_vars_file) {
  tbls = lapply(vars_files,function(x) {
    tbl = data.table::fread(x$out_vars_pos)
    tbl[,subtype:=x$subtype]
    tbl[,gene:=x$gene]
  })
  tbl = data.table::rbindlist(tbls)
  tbl_col = quantcut.ordered.color(tbl$freq,as.rgb = T)
  ordered.color.legend.report(tbl_col)
  tbl = cbind(tbl,tbl_col[,.(red,green,blue)])
  write.table(tbl,out_vars_file,row.names = F,sep="\t")
  report$add.descr(sprintf("We have created the file '%s' to color-code the variability
                           on the 3D structure. You need to start Pymol from the directory
                           that holds that file; source 'rsv_surveilance_utils.py' with the
                           'run' command; and execute 'remap_and_paint()' function. See
                           comments in the Python source for the arguments.",
                           out_vars_file))
  tbl
}

rsv_ali_to_index <- function(ali_file,id_from=NULL,id_to,index_file) {
  ## this index is used to remap residue numbers in the PDB file - 
  ## this works for RSV F because there are no insertions in the
  ## implicit PDB sequence relative to the RADAR reference
  ## Normally, sequence IDs in the alignment (except for id_to) will
  ## be PDB chain IDs.
  ali = readAAMultipleAlignment(filepath = ali_file,format="fasta")
  if(is.null(id_from)) {
    id_from = setdiff(rownames(ali),id_to)
  }
  ali.d = as.data.table(t(as.matrix(ali)))
  res = NULL
  for(id_from_one in id_from) {
    ali.d.one = copy(ali.d[,c(id_from_one,id_to),with=F])
    setnames(ali.d.one,c("from","to"))
    ali.d.one[,`:=`(id_from=id_from_one,
                    id_to=id_to,
                    from_ind=(.I-cumsum(from=="-")),
                    to_ind=(.I-cumsum(to=="-")))]
    res = rbind(res,ali.d.one)
  }
  write.table(res,index_file,row.names = F,sep="\t")
  res
}

rsv_tree_report <- function(ali_file) {
  ali = readAAMultipleAlignment(filepath = ali_file,format="fasta")
  #make.global()
}

new_ordination.task <- function(main.meta.var,norm.method,label=NULL,size=NULL,distance=NULL) {
  distance.arg = distance
  ord.method = "NMDS"
  distance.0="bray"
  within(mgsat.16s.task.template$test.counts.task$ordination.task, {
    distance=if(is.null(distance.arg)) distance.0 else distance.arg
    ord.tasks = list(
      list(
        ordinate.task=list(
          method=ord.method
          ##other arguments to phyloseq:::ordinate
        ),
        plot.task=list(
          type="samples",
          color=main.meta.var,
          label = label,
          size = size,
          axis.scale = c(1,1,1)
          ##other arguments to phyloseq:::plot_ordination
        )
      ),
      list(
        ordinate.task=list(
          method=ord.method,
          formula=main.meta.var
          ##other arguments to phyloseq:::ordinate
        ),
        plot.task=list(
          type="samples",
          color=main.meta.var,
          label = label,
          size = size,
          axis.scale = c(1,1,1)
          ##other arguments to phyloseq:::plot_ordination
        )
      )          
    )
  })            
}


load_rsvdb_seq <- function(seq.file,meta.file,is.ali=T,dedup=T,clean.seq=T,sanitize.seq.id=T,
                           us_states_file=NULL,ali.file=NULL,
                           msa.method = "MAFFT",
                           msa.args = list(param="--maxiterate 1000 --localpair --thread -1"),
                           preclust.ident=0.99) {
  library(Biostrings)
  library(tidyr)
  stopifnot(file.exists(seq.file))
  if(is.ali) {
    s = as(readAAMultipleAlignment(seq.file,format="fasta"),"XStringSet")
  }
  else {
    s = Biostrings::readAAStringSet(seq.file)
  }
  if(clean.seq) {
    s = remove.end.stop.seq(s,chop.all = T)
  }
  nm = names(s)
  seq.id = id.from.seq.names(nm)
  seq.id.orig = seq.id
  if(sanitize.seq.id) {
    seq.id = sanitize.seq.id.for.phangorn(seq.id)
  }
  seq.id.dt = data.table(SeqIDSanitized=seq.id,SeqID=seq.id.orig,SeqName=nm)
  meta = data.table::fread(meta.file)
  meta = meta[seq.id.dt,on="SeqID",nomatch=0]
  stopifnot(all(meta$SeqID==seq.id.dt$SeqID))
  meta[,`:=`(SeqIDOrig=SeqID,SeqID=SeqIDSanitized,SeqIDSanitized=NULL)]
  meta[,rn:=SeqID]
  meta[,State:=meta[,stringr::str_split_fixed(OUTSMART_ID,"_",2)[,1]]]
  meta[,Age_month_qua:=quantcut.ordered(Age_month)]
  meta[,Hospital_Stay_Long:=ifelse(grepl("^>24H",Hospital_Stay),">24H","<24H")]
  if(!is.null(us_states_file)) {
    us_states = data.table::fread(us_states_file)
    ## https://www.census.gov/geo/reference/gtc/gtc_census_divreg.html
    ## Puerto Rico is not part of any census region, will get NA
    meta = us_states[,.(State=STUSPS,RegionID=REGION,Region=REGION_NAME,StateName=NAME)][meta,on="State"]
    ## Reset NA region to full StateName, but keep RegioID as 9 as it is in the us_states shape file
    ## for whatever reason
    meta[is.na(Region),Region:=StateName]
    meta[,Region:=factor(Region,levels = c("West","Midwest","South","Northeast","Puerto Rico"))]
  }
  stopifnot(all(names(s)==meta$SeqName))
  meta$n_nongenic = nongenic.frequency.seq(s,out.format = "rowSums")
  names(s) = meta$SeqID
  seq = list(seq=s,attr=meta)
  seq = init_seq_for_dedup(seq)
  if(dedup) {
    seq = dedup_seq(seq,attr_aggr_fun = rsv_aggr_attr_exemplar)
    aggregation.counts.report(seq)
  }
  if(file.dep.updated(seq.file,ali.file)) {
    ali = do.call(mgsat.msa,
                  c(list(seq$seq,msa.method=msa.method),msa.args))
    seq$seq = as(ali,"XStringSet")
    dir.create(dirname(ali.file))
    writeXStringSet(seq$seq,ali.file,format = "fasta")
  }
  else {
    seq$seq = as(readAAMultipleAlignment(ali.file,format="fasta"),"XStringSet")
  }
  stopifnot(all(names(seq$seq)==seq$attr$SeqName))
  seq = aggregate.seq.by.sim(seq,preclust.ident = preclust.ident)
  stopifnot(all(names(seq$seq)==seq$attr$SeqName))
  return (seq)
}

aggregate.seq.by.sim <- function(seq,
                                 n.seq.max = 200,
                                 preclust.ident = 0.99) {
  
  ali = Biostrings::AAMultipleAlignment(seq$seq)
  
  cluster.dist.type = "DECIPHER.sim"
  ali.der = rsv.ali.deriv(ali,type=cluster.dist.type)
  ## if there are too many sequences, cluster them and mask everything but selected exemplars
  n.seq = length(labels(ali.der$d))
  if(n.seq > n.seq.max) {
    cl.res = cluster.with.dist.cutoff(d=ali.der$d,
                                      h=ident.to.seq.dist(preclust.ident,dist.transform = ali.der$dist.transform))
    exemplars = pick.dist.group.exemplar(d=ali.der$d,group=cl.res$group,seq.attr=seq$attr)
    exemplar.report(exemplars)
    seq = aggregate_seq(seq=seq,exemplars = exemplars$dind,attr_aggr_fun = rsv_aggr_attr_exemplar)
    report$add.descr(sprintf("**Pre-clustered (by complete linkage) with identity cutoff %s from %s sequences into %s exemplars**",
                             preclust.ident,n.seq,nrow(seq$attr)))
    aggregation.counts.report(seq)
  }
  seq
}

aggregation.counts.report <- function(seq) {
  report$add.table(seq$attr[,.(Count=sum(Count),N_uniq=.N),by=Subtype],caption = "Counts of nonredundant sequences")
  report$add.table(seq$attr_all[,.(Count=sum(Count),N_uniq=.N),by=Subtype],caption = "Counts of all sequences")
}

plot.rsv.geo <- function(meta.seqs,geo_db_dir) {
  library(tmap)
  attr_all = meta.seqs$attr_all
  #gr_summ = rsv_aggr_attr(attr_all,by="RegionID")
  gr_summ = rsv_aggr_attr(attr_all,by=c("RegionID","State"))
  report$add.table(gr_summ,caption = "Summary table used for geomaps")
  #us_regions = readRDS(file.path(geo_db_dir,"us_regions.rds"))
  us_states = readRDS(file.path(geo_db_dir,"us_states.rds"))
  #us_regions_cnt = tmaptools::append_data(us_states_cont,gr_summ, key.shp = "REGION", key.data = "RegionID")
  us_states_summ = tmaptools::append_data(us_states,gr_summ, key.shp = "STUSPS", key.data = "State")
  us_states_cont = us_states_summ[us_states_summ$RegionID %in% c("1","2","3","4") & us_states_summ$STUSPS != "AK",]
  us_regions = maptools::unionSpatialPolygons(us_states, IDs=us_states$REGION)
  
  us_states_cont$Label=sprintf("%s\n(%s)",us_states_cont$STUSPS,us_states_cont$Count)
  
  #make.global()
  
  for(var_opt in list(list(col="Count"),
                      list(col="Hospital_Stay_Long_Ratio",
                           style="fixed",breaks=seq(0.0,1.0,by=0.2),
                           interval.closure="right"),
                      list(col="Age_group_median"))) {
    ##is.master to tm_shape can be used to constraint the bounding box to filtered regions,
    ##without the need to subset the states?
    pl = tm_shape(us_states_cont, projection=2163) + 
      do.call(tm_polygons,c(
        list(border.col = "grey50", 
             border.alpha = .5, 
             title = "", 
             showNA = F),
        var_opt)) +
      # tm_polygons(var_opt$var_name, 
      #             border.col = "grey50", 
      #             border.alpha = .5, 
      #             title = "", 
      #             showNA = F,
      #             style=var_opt$style,
      #             breaks=var_opt$breaks,
      #             interval_closure=var_opt$interval_closure) +
      # 
      tm_text("Label", col = "grey30", root=3) +
      tm_shape(us_regions) +
      tm_borders(lwd=1, col = "black", alpha = .5) +
      tm_credits("Shapes @ Unites States Census Bureau", position = c("left", "bottom")) +
      tm_layout(scale = 1.5,fontface="bold",legend.text.size = 1)
    
    report$add(pl,
               caption = sprintf("Choropleth map for %s",var_opt$col))
    
  }
}

plot.rsv.tree.make.tip.insets <- function(tr.df,attr_all,inset.type=c("bar_stacked","pie")) {
  inset.type=inset.type[1]
  tr.df = tr.df[,.(SeqID,node)]
  attr_all = attr_all[tr.df,on=c(SeqID_uniq="SeqID")]
  attr_all[,Region:=factor(Region)]
  pal = generate.colors.mgsat(attr_all[,Region],value = "palette")
  lapply(split(attr_all,by="node"),function(dat) {
    if(inset.type=="pie") {
      pl.ins = ggplot(dat,aes(x=Count,fill=Region)) +
        geom_bar(width=1,alpha=0.5)+coord_polar(theta="y")
    }
    else if(inset.type=="bar_stacked") {
      pl.ins  = ggplot(dat,aes(x=factor(1),fill=Region)) + 
        geom_bar(position="stack",color="black",alpha=0.5) + coord_flip()
    }
    pl.ins = pl.ins + scale_fill_manual(values=pal) + scale_colour_manual(values=pal) + ggtree::theme_inset()
  }
  )
}

srv.region.cum_cnt <- function(attr_all) {
  seq_cnt = attr_all[,.N,by=.(SeqID_uniq,Region)][order(SeqID_uniq,Region)][,.(Region,end=cumsum(N),ind_reg=seq(.N)),by=SeqID_uniq]
  seq_cnt[,iid:=.I]
  seq_cnt[,iid_p:=ifelse(ind_reg>1,.I-1,.I)]
  seq_cnt[,start:=ifelse(ind_reg>1,end[iid_p]+1,1)]
  seq_cnt[,.(SeqID_uniq,Region,start,end)]
}

mgsat_theme_transparent <- function (...) 
{
  theme_bw() + theme(panel.background = element_rect(fill = "transparent", 
                                        colour = NA), 
        plot.background = element_rect(fill = "transparent", 
                                       colour = NA), 
        legend.key = element_rect(fill = "transparent", 
                                  colour = NA), 
        legend.background = element_rect(fill = "transparent", 
                                         colour = NA), ...)
}

mgsat_theme_transparent <- function (base_size = 11, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% theme(panel.background = element_blank(), 
                     plot.background = element_blank(), 
                     legend.key = element_blank(), 
                     legend.background = element_blank(), 
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     panel.border=element_blank(),
                     axis.line.x=element_blank(),
                     axis.line.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.ticks.y=element_blank(),
                     complete = TRUE)
}

mgsat_theme_transparent_tree_with_scale <- function(...) {
  mgsat_theme_transparent(...) %+replace% theme(axis.line.x=element_line(colour="black"),
                                       axis.ticks.x=element_line(colour="black"),
                                       axis.text.x=element_text(colour="black",size=rel(0.8),margin = margin(3)),
                                       complete=TRUE)
}

mgsat_facet_plot <- function (p, panel, data, geom, mapping = NULL, ...) 
{
  p <- add_panel(p, panel)
  df <- p %+>% data
  p + geom(data = df, mapping = mapping, ...)
}

#at 97% clustering:
#region_stacked_width=600
#region_stacked_height=800
plot.rsv.tree <- function(tree.attr,meta.seqs,gene,subtype,count_scale_limits=NULL,
                          region_stacked_width=800,
                          region_stacked_height=1000) {
  caption = sprintf("Tree for gene %s subtype %s",gene,subtype)
  report.section = report$add.header(
    caption,
    section.action = "push")
  
  library(ggtree)
  
  tree = tree.attr$tree.dist$tree
  
  if(subtype=="AB") {
    tree = root(tree,outgroup=which(grepl("A$",tree$tip.label)),r=T)  
  }
  seq.attr = tree.attr$seq.attr
  setkey(seq.attr,rn)
  
  top_SeqID = meta.seqs$attr[Count==max(Count),SeqID]
  top_attr = meta.seqs$attr_all[SeqID_uniq %in% top_SeqID]
  top_attr = top_attr[,.(Count=sum(Count)),by=.(SeqID_uniq,State)][order(Count,decreasing = T)]
  
  report$add.table(top_attr,caption = "Distribution of the most abundant sequence(s) across states")
  
  report$add.table(top_attr[,.(CountStates=.N),by=SeqID_uniq],caption = "Number of states with the most abundant sequence(s)")
  
  priv.add.tree.plot.report <- function(branch.length="edge.length",count_scale_limits_p=count_scale_limits) {
    
    ##geom_point2 gives an error; the only method to plot a specific subset that worked is below:
    #tree_layout="circular"
    tree_layout="rectangular"
    pl.df = data.table(fortify(tree,layout = tree_layout,branch.length = branch.length),
                       keep.rownames = T)
    pl.df[,rn:=label]
    pl.df = pl.df[!is.na(rn),]
    setkey(pl.df,rn)
    #tr.df = seq.attr[tree$tip.label,nomatch=0]
    tr.df = seq.attr[pl.df,nomatch=0]
    ## first column of the user data matrix must match the tip labels
    setcolorder(tr.df,c("rn",colnames(tr.df)[colnames(tr.df)!="rn"]))
    #rownames(tr.df) = tr.df$rn
    
    if(branch.length=="edge.length") {
      descr = "phylogram (with branch lengths)"
    }
    else {
      descr = "cladogram (without branch lengths)"
    }
    pal = generate.colors.mgsat(meta.seqs$attr_all[,Region],value = "palette")
    ##NJ tree aims to minimize total tree length (sum of all branches). It is a greedy algorithm
    ##that at every step maximizes the reduction in tree length. Thus, we can see consensus and
    ##nearconsensus in different branches.
    #size=ident.mean
    pl_base = ggtree(tree,layout=tree_layout,branch.length = branch.length,size=0.8) +
      mgsat_theme_transparent_tree_with_scale(base_size=22) + 
      theme(legend.position="right",axis.line.x=element_line(colour="black")) +
      scale_x_continuous(expand=c(0,0))
    reg_cnt = as.data.frame(dcast(meta.seqs$attr_all,SeqID_uniq ~ Region,value.var = "Count",fun.aggregate = function(x) sum(x,na.rm = T),fill=NA))
    rownames(reg_cnt) = reg_cnt$SeqID_uniq
    reg_cnt$SeqID_uniq = NULL
    skip_some = T
    if(!skip_some) {
      pl = gheatmap(pl_base,reg_cnt,color="black",font.size = 6) + scale_fill_gradientn(colours=rainbow(10,s=0.3,v=0.7,start=2/6,end=1),na.value = "white")
      report$add(pl,
                 caption = sprintf("Phylogenetic tree 
                                 in circular %s layout.
                                 Heatmap represents counts of identical sequences across US Census Regions and Puerto Rico", descr),
                 width=1000,
                 height=800)
      pl = pl + coord_polar(theta='y')
      report$add(pl,
                 caption = "Same as above with polar coordinates")
      pl.ins = plot.rsv.tree.make.tip.insets(tr.df=tr.df,attr_all=meta.seqs$attr_all,inset.type = "pie")
      pl = ggtree::inset(pl_base,pl.ins,height = 0.05)
      report$add(pl,
                 caption = "Alternative represention of sequence counts with pie charts")
    }
    reg_cnt_cum = srv.region.cum_cnt(meta.seqs$attr_all)
    ##ATTENTION: facet_plot internally uses `%+>%` operator, which calls merge on
    ##the user data frame and a data frame built from tree tip labels and y positions,
    ##and uses first column of the user data frame as key (column name does not matter).
    ##Below, we make this explicit by creating the first column named id. Note that in our case,
    ##the internal facet_plot merge is a many-to-one operation (several regions per each unique id)
    reg_cnt_cum = data.table(id=reg_cnt_cum$SeqID_uniq,reg_cnt_cum)
    ##this is not needed - see comment above
    #reg_cnt_cum_j = tr.df[,.(SeqID_uniq=label,y)][reg_cnt_cum,on="SeqID_uniq"]
    #stopifnot(all(reg_cnt_cum_j$SeqID_uniq==reg_cnt_cum$SeqID_uniq))
    #stopifnot(!reg_cnt_cum_j[,any(is.na(y))])
    #reg_cnt_cum = reg_cnt_cum_j
    panel_name = "Number of Isolates"
    pl = facet_plot(pl_base,data=reg_cnt_cum,panel=panel_name, geom_segment, 
                    aes(x=start-1, xend=end+1e-4, y=y, yend=y,color=Region),size=rel(2))
    pl = pl + scale_fill_manual(values=pal) + scale_colour_manual(values=pal)
    report$add(pl,
               caption = "Alternative represention of sequence counts with stacked bars",
               width=region_stacked_width,
               height=region_stacked_height)
    pl_proto_multiyear = facet_plot(pl,data=reg_cnt_cum,panel="Number of Isolates 2018", geom_segment, 
                                    aes(x=start-1, xend=end+1e-4, y=y, yend=y,color=Region),size=rel(2))
    report$add(pl_proto_multiyear,
               caption = "Crude prototype for showing multi-season distribution of sequences",
               width=region_stacked_width,
               height=region_stacked_height)    
    ##once facet_plot (add_panel) has been called on the plot, we can use the operator to
    ##extract tree coords into our data frame
    #reg_cnt_cum_tr = as.data.table(pl  %+>% reg_cnt_cum)
    #lab.df = reg_cnt_cum_tr[,.(end=max(end),y=y[1]),by=label]
    reg_cnt_cum_maxend = reg_cnt_cum[,.(max_end=max(end)),by=id]
    reg_cnt_cum = reg_cnt_cum[reg_cnt_cum_maxend,on="id"]
    ## Note that facet_plot is implemented with facet_grid on a hidden variable that has 'panel' values
    ## as levels, so adding another geom to the existing is as simple as calling facet_plot again with
    ## the same value of 'panel'.
    ## hjust="inward" prevents clipping of labels; other alternative is here with ggplot_gtable
    ## https://stackoverflow.com/questions/17241182/how-to-make-geom-text-plot-within-the-canvass-bounds
    pl_lab = facet_plot(pl,data=reg_cnt_cum,panel=panel_name, geom_text, 
                                    aes(x=max_end, y=y,label=label),size=rel(2),hjust="inward")
    report$add(pl_lab,
               caption = "Large alternative represention of sequence counts with stacked bars and tip labels, for reference",
               width=region_stacked_width,
               height=region_stacked_height*2)    
    # geom_tiplab() +
    # geom_point(data=pl.df,aes(size=ident.mean,fill=descr),color="black",shape=21,alpha=1) +
    # geom_text(data=pl.df,aes(label=paste(descr,rn),color=descr),nudge_x=0,nudge_y=0,alpha=1) +
    if(!skip_some) {
      pl = pl_base %<+% tr.df + 
        geom_tippoint(aes(color=Hospital_Stay_Long_Ratio,size=Count),alpha=0.7) +
        # this scales area
        scale_size(range = c(4, 18),limits = count_scale_limits_p) + 
        scale_color_gradient2(midpoint = mean(tr.df$Hospital_Stay_Long_Ratio,na.rm=T),na.value = "lightgrey") +
        coord_polar(theta='y')
      report$add(pl,
                 caption = "Tree colored by percentage of long (>24H) hospital stay")
      
      pal_age_group = generate.colors.mgsat(tr.df[,Age_group_median],value = "palette")
      
      pl = pl_base %<+% tr.df + 
        geom_tippoint(aes(color=Age_group_median,size=Count),alpha=0.7) +
        # this scales area
        scale_size(range = c(4, 18),limits = count_scale_limits_p) +
        coord_polar(theta='y') +
        scale_fill_manual(values=pal_age_group) + scale_colour_manual(values=pal_age_group)
      report$add(pl,
                 caption = "Tree colored by the mean age")
    }
    tr.df
  }
  for(branch.length in c("edge.length")) {
    #for(branch.length in c("edge.length","none")) {
    pl.df = priv.add.tree.plot.report(branch.length,count_scale_limits_p = count_scale_limits)
  }
  report$add.table(pl.df,caption = "Data used to plot the tree")
  report$pop.section()
}

plot.rsv.ordinations <- function(tree.attr) {
  seq.attr = tree.attr$seq.attr
  seq.attr[,SampleID:=SeqID]
  setkey(seq.attr,rn)
  dist.mat = as.dist(tree.attr$d.mat[seq.attr$rn,seq.attr$rn])
  m_a = list(attr=seq.attr)
  ordination.task = new_ordination.task("Age_group_median",label="SeqID",distance=dist.mat)
  do.call(ordination.report,
          c(list(m_a=m_a,res=list(),descr="Colored by group of the median age"),
            ordination.task
          )
  )  
  ordination.task = new_ordination.task("Hospital_Stay_Long_Ratio",label="SeqID",distance=dist.mat)
  do.call(ordination.report,
          c(list(m_a=m_a,res=list(),descr="Colored by ratio of long hospital stay"),
            ordination.task
          )
  )  
  
}

rsv.download.geo <- function(db_dir) {
  library(data.table)
  us_regions = tigris::regions(resolution = '20m')
  us_states = tigris::states()
  saveRDS(us_states,file.path(db_dir,"us_states.rds"))
  saveRDS(us_regions,file.path(db_dir,"us_regions.rds"))
  st = data.table(us_states@data)
  reg = data.table(us_regions@data)
  st = reg[,.(REGION_NAME=NAME,REGION=REGIONCE)][st,on="REGION"]
  write.table(st,file.path(db_dir,"us_states.csv"),sep=",",row.names = F)
}

prepare_tree_input <- function(ali,ali_file) {
  seq = as(ali,"XStringSet")
  ## RaxML complaints when not all letters from the substitution
  ## model are found in the sequence. Although it proceeds, this
  ## post suggests that the results might be faulty:
  ## https://groups.google.com/forum/#!topic/raxml/yHvXAKdk7OA
  ## Therefore, we append columns of missing characters
  ## to the alignment. Because they are identical in all sequences,
  ## this should not change the tree
  fr = seq_letter_frequency(seq)
  miss = fr$not_in_seq
  if(length(miss)>0) {
    ## xscat acts like paste0 and recycles the arguments,
    ## so this will append miss string at the end of every seq
    seq_names = names(seq)
    seq = Biostrings::xscat(seq,paste0(miss,collapse = ""))
    names(seq) = seq_names
  }
  dir.create(dirname(ali_file),recursive = T)
  Biostrings::writeXStringSet(seq,ali_file)
}

make_raxml_tree <- function(tree_dir,ali_file,id_task,threads=4) {
  sfx = id_task
  #model = "PROTCATAUTO"
  #model = "PROTGAMMALG4X"
  model = "PROTGAMMAAUTO"
  #model = "PROTGAMMAHIVWF"
  cmd = sprintf("cd %s && raxmlHPC-PTHREADS-SSE3 -T %s -f a -m %s -p 12345 -x 12345 -# autoMRE -s %s -n %s",
                tree_dir,threads,model,normalizePath(ali_file),sfx)
  message(sprintf("Building ML tree with command: `%s`",cmd))
  system(cmd)
}

rsv_process <- function() {
  report.section = report$get.section()
  report$add.header("Iterating over genes and subtypes")
  report$push.section(report.section)
  do.pdb.only = F
  do.plots = T
  do.plots.tree = T
  do.plots.ord = F
  do.geo = F
  do.tree_build_ml = F
  
  cons_freq_type = "CONS"
  
  tree.method = "nj" #"load" #"nj"
  tree.dist.type = "DECIPHER.sim"
  
  outlier.k=11
  msa.method = "MAFFT"
  msa.args = list(param=sprintf("--maxiterate 1000 --localpair --thread %s",options("mc.cores")))  
  
  cluster.dist.type = "DECIPHER.sim"
  preclust.ident=0.97
  
  ## A kludge to set point size on the tree plot to uniform scale across subtypes
  ## when scaling points to the count of each unique sequence. Set it first to NULL,
  ## run for all subtypes, look up automatically computed ranges and set her to the union
  ## across subtypes.
  ## Note: it appears that the limits should be (min,max) count; if it is narrower,
  ## we get all dot-sized points
  
  count_scale_limits = c(1,90)
  
  season = "Season_2017"
  db_dir = file.path("../../../data",season,"database")
  data_dir = file.path("../../../data",season,"OSMT_2017-11-20")
  
  geo_db_dir = file.path(db_dir,"geo")
  #rsv.download.geo(geo_db_dir)
  
  meta_file = file.path(data_dir,"meta.csv")
  nuc_ali_dir = file.path(data_dir,"nuc_ali")
  nuc_files_patt = file.path(data_dir,"*.\\.fasta$")
  prot_dir = file.path(data_dir,"protein")
  prot_ali_dir = file.path(data_dir,"prot_ali")
  #rsv_nuc_to_prot(nuc_files_patt = nuc_files_patt,nuc_ali_dir = nuc_ali_dir)
  pdb_out_dir = file.path(data_dir,"pdb_out")
  pdb_inp_dir = file.path(data_dir,"PDB")
  radar_out_dir = file.path(data_dir,"RADAR_Output")
  tree_dir = file.path(data_dir,"tree")
  
  pdb_seqs = file.path(pdb_inp_dir,c("4jhw_f.fasta","3rrr_f.fasta"))
  pdb_vars = file.path(pdb_out_dir,"RSV_F_vars.txt")
  ## Process note: half of BG sequences are truncated at 3', started at the
  ## last base of AAT (or AAC) codon (10 bases total - see the aligned nuc sequences). 
  ## To remove noise from the tree, we truncate
  ## all BG sequences in the protein MSA to the last intact AA.
  ## The truncation coordinates are communicated through 
  tasks = list(
    list(gene="F",subtype="A",
         id_task = "Osmt2017_RSV_AF",
         ref_radar_id="ref_FA_KX858757",
         ref_prot_seq=file.path(radar_out_dir,"ref_FA_KX858757.prot.fasta"),
         radar_var=file.path(radar_out_dir,"OS2_RSVAF_VR_CDS.txt"),
         prot = file.path(prot_dir,"Osmt2017_RSV_AF.prot.fasta"),
         prot_ali = file.path(prot_ali_dir,"Osmt2017_RSV_AF.prot.ali.fasta"),
         prot_ali_range = c(1,574), #cuts ambigous J; TODO: replace all ambiguous with X
         prot_ali_tree = file.path(tree_dir,"Osmt2017_RSV_AF.prot.ali.fasta") 
    ),
    list(gene="F",subtype="B",
         id_task = "Osmt2017_RSV_BF",
         ref_radar_id="ref_FB_KX858756",
         ref_prot_seq=file.path(radar_out_dir,"ref_FB_KX858756.prot.fasta"),
         radar_var=file.path(radar_out_dir,"OS2_RSVBF_VR_CDS.txt"),
         prot = file.path(prot_dir,"Osmt2017_RSV_BF.prot.fasta"),
         prot_ali = file.path(prot_ali_dir,"Osmt2017_RSV_BF.prot.ali.fasta"),
         prot_ali_tree = file.path(tree_dir,"Osmt2017_RSV_BF.prot.ali.fasta") 
    ),
    list(gene="F",subtype="AB",
         id_task = "Osmt2017_RSV_F",
         ref_radar_id=NULL,
         ref_prot_seq=NULL,
         radar_var=NULL,
         prot = file.path(prot_dir,"Osmt2017_RSV_F.prot.fasta"),
         prot_ali = file.path(prot_ali_dir,"Osmt2017_RSV_F.prot.ali.fasta"),
         prot_ali_range = c(1,574), #cuts ambigous J; TODO: replace all ambiguous with X
         prot_ali_tree = file.path(tree_dir,"Osmt2017_RSV_F.prot.ali.fasta")         
    ),
    list(gene="G",subtype="A",
         id_task = "Osmt2017_RSV_AG",
         prot = file.path(prot_dir,"Osmt2017_RSV_AG.prot.fasta"),
         prot_ali = file.path(prot_ali_dir,"Osmt2017_RSV_AG.prot.ali.fasta"),
         prot_ali_tree = file.path(tree_dir,"Osmt2017_RSV_AG.prot.ali.fasta")         
    ),
    list(gene="G",subtype="B",
         id_task = "Osmt2017_RSV_BG",
         prot = file.path(prot_dir,"Osmt2017_RSV_BG.prot.fasta"),
         prot_ali = file.path(prot_ali_dir,"Osmt2017_RSV_BG.prot.ali.fasta"),
         prot_ali_range = c(1,107),
         prot_ali_tree = file.path(tree_dir,"Osmt2017_RSV_BG.prot.ali.fasta")
    ),
    list(gene="G",subtype="AB",
         id_task = "Osmt2017_RSV_G",
         prot = file.path(prot_dir,"Osmt2017_RSV_G.prot.fasta"),
         prot_ali = file.path(prot_ali_dir,"Osmt2017_RSV_G.prot.ali.fasta"),
         #prot_ali_range = c(1,122),
         prot_ali_tree = file.path(tree_dir,"Osmt2017_RSV_G.prot.ali.fasta")
    )
  )
  tasks_upd = list()
  for(task in tasks) {
    report$add.header(sprintf("Analysis for gene %s subtype %s",task$gene,task$subtype))
    report$push.section(report.section)
    tasks_pdb = list()
    if(!is.null(task$radar_var)) {
      for(pdb_seq in pdb_seqs) {
        task_pdb = list()
        task_pdb$pdb_radar_ref_ali = rsv_pdb_radar_ref_ali(pdb_seq = pdb_seq,
                                                           ref_seq = task$ref_prot_seq,
                                                           pdb_out_dir = pdb_out_dir)
        task_pdb$pdb_renumber_index = replace_file_ext(task_pdb$pdb_radar_ref_ali,".txt")
        rsv_ali_to_index(ali_file=task_pdb$pdb_radar_ref_ali,
                         id_to=task$ref_radar_id,
                         index_file=task_pdb$pdb_renumber_index)
        tasks_pdb = c(tasks_pdb,list(task_pdb))
      }
      task$tasks_pdb = tasks_pdb
      task$out_vars_pos = rsv_pdb_selections(radar_var=task$radar_var,pdb_out_dir=pdb_out_dir,freq_type = cons_freq_type)
    }
    if(!do.pdb.only) {
      meta.seqs = load_rsvdb_seq(seq.file=task$prot,
                                 meta.file=meta_file,
                                 is.ali=F,
                                 dedup=T,
                                 clean.seq=T,
                                 sanitize.seq.id=T,
                                 us_states_file=file.path(geo_db_dir,"us_states.csv"),
                                 ali.file=task$prot_ali,
                                 msa.method = msa.method,
                                 msa.args = msa.args,
                                 preclust.ident=preclust.ident
      )
      #make.global()
      ali = Biostrings::AAMultipleAlignment(meta.seqs$seq)
      alignment.report(ali,caption = sprintf("Multiple alignment of all %s sequences",nrow(ali)))
      if(!is.null(task$prot_ali_range)) {
        ali = Biostrings::AAMultipleAlignment(as(ali,"XStringSet"),
                                              start=task$prot_ali_range[1],
                                              end=task$prot_ali_range[2])
        alignment.report(ali,
                         caption = sprintf("Multiple alignment of all %s sequences **after** trimming the 
                                           alignment to a range (%s)",
                                           nrow(ali),
                                           paste(task$prot_ali_range,collapse = ",")))
        rsv.sequence.qc.report(as(ali,"XStringSet"),meta.seqs$attr)
        ##TODO: dedup here again since trimming could have created identical sequences
      }
      if(do.tree_build_ml) {
        prepare_tree_input(ali,task$prot_ali_tree)
        make_raxml_tree(tree_dir=tree_dir,ali_file=task$prot_ali_tree,id_task=task$id_task,threads=options("mc.cores"))
      }
      
      if(tree.method=="load") {
        ## ggtree::read.raxml can only read boostraps on the edges (RAxML_bipartitionsBranchLabels.* file),
        ## while FigTree can only read them on the nodes (RAxML_bipartitions.* file) 
        task$prot_tree = file.path(tree_dir,sprintf("RAxML_bipartitionsBranchLabels.%s",task$id_task))
      }
      
      ali.der = rsv.ali.deriv(ali,type=cluster.dist.type)
      
      report$add.header("Checking for the outliers",sub=T,report.section = report.section)
      
      outlier.ordinate.method = "NMDS"
      if(rsv.vac.debug.mode) {
        outlier.ordinate.method = "MDS"
      }
      outl.res = mgsat.find.outliers.ordinate(ali.der$d,
                                              k = outlier.k,
                                              alpha = 0.01,
                                              pval.adjust="BH",
                                              lof.on.ordinate = F,
                                              ordinate.args=list(method=outlier.ordinate.method,k=3))
      #keep or remove the outliers
      outl_names = names(outl.res$ind.outlier)
      ##because we already built the ML tree above, we should not remove rows from the alignment
      ##afterwards. TODO: move ML tree construction into rsv.ali.to.tree.attr, so that row and column
      ##masking is respected
      if(! tree.method=="load") {
        rowmask(ali,append="union") = IRanges(start=(rownames(ali) %in% outl_names))
      }
      if(length(outl_names)>0) {
        report$add.table(meta.seqs$attr[SeqID %in% outl_names],
                         caption = "Metadata for detected outlier sequences")
      }
      report$pop.section()
      
      meta.seqs$attr = filter.seq.attr.by.ali(ali,meta.seqs$attr,drop.rows=T)
      
      report$add.header("Removing any all-gap columns")
      ## make sure that we will not drop any previously masked columns 
      ## that are not all gaps
      colmask(ali) = NULL
      ##mask.gaps.ali respects set row masks
      ali = mask.gaps.ali(ali, min.fraction=1.0, min.block.width=1)
      ## get the masked rows and columns out of the way for good.
      ali = masked.copy.ali(ali,drop.rows = T,drop.columns = T,deep = T)
      meta.seqs$seq = ali
      
      report$add.table(meta.seqs$attr,
                       caption = "Attributes of the unique sequence set")
      report$add.table(meta.seqs$attr[,.(SeqID,Count,SeqID_all)],
                       caption = "Selected representative IDs of the identical
                       sets of sequences (including singleton clusters)")
      
      if(do.plots) {
        tree.attr = rsv.ali.to.tree.attr(ali=ali,
                                         seq.attr=meta.seqs$attr,
                                         use.tree.dist = F,
                                         make.tree = T,
                                         mask.gaps = T,
                                         root.tree = T,
                                         tree.method = tree.method,
                                         ##ATTENTION:
                                         ##We knock off all gaps because it this dataset there is 
                                         ##a single genotype, and gaps are artefacts.
                                         ##Change this in future seasons if gaps become meaningful.
                                         mask.gaps.min.fraction=0.00001,
                                         mask.gaps.min.block.width=1,
                                         tree.file=task$prot_tree,
                                         tree.dist.type=tree.dist.type)
        if(do.plots.tree) {
          plot.rsv.tree(tree.attr,meta.seqs=meta.seqs,
                        gene=task$gene,task$subtype,
                        count_scale_limits = count_scale_limits)
        }
        if(do.plots.ord) {
          plot.rsv.ordinations(tree.attr)
        }
        if(do.geo) {
          report$add.header(sprintf("Geomaps for gene %s subtype %s",task$gene,task$subtype))
          report$push.section(report.section)
          plot.rsv.geo(meta.seqs=meta.seqs,geo_db_dir=geo_db_dir)
          report$pop.section()
        }
        
      }
    }
    tasks_upd = c(tasks_upd,list(task))
    report$pop.section()
  }
  vars_files = tasks_upd[sapply(tasks_upd,function(x) !is.null(x$out_vars_pos))]
  if(length(vars_files)>0) {
    rsv_make_vars_colors(vars_files=vars_files,
                         out_vars_file=pdb_vars)
  }
  report$pop.section()
}

rsv.install.packages <- function() {
  ##'EBIImage', which is a dependency of ggtree, needs fftw3 library (e.g. `brew install fftw`)
  bio.pkgs = c("Biostrings","ggtree","DECIPHER")
  vanilla.pkgs = c("seqinr","picante","ape","curl","phangorn")
  install.packages(vanilla.pkgs)
  source("http://bioconductor.org/biocLite.R")
  biocLite(bio.pkgs) 
  ## rMSA depends on proxy
  devtools::install_version("proxy") 
  #conda install mafft
  devtools::install_github("mhahsler/rMSA")
  devtools::install_github("tidyverse/tibble") # plotly fails with CRAN version with missing as_tibble()
}

rsv.install.geo.packages <- function() {
  #tigris needs gdal libraries in the OS
  #brew install gdal or conda install gdal
  #conda install jq r-udunits2
  install.packages(c("tigris"))
  #brew install homebrew/versions/v8-315
  devtools::install_github("cran/acs")
  devtools::install_github("mtennekes/tmaptools")
  devtools::install_github("mtennekes/tmap")
}

threads = 4

## number of cores to use on multicore machines
options(mc.cores=threads)
options(boot.ncpus=threads)
## parallel backend
options(boot.parallel="snow")


## location of MGSAT code
MGSAT_SRC = "~/work/mgsat"

source(paste(MGSAT_SRC,"dependencies.r",sep="/"),local=T)

## Uncomment next line to install packages needed by MGSAT (!!!comment it out
## in all subsequent runs once the packages have been installed!!!).
## Note: you should also pre-install Pandoc program from http://johnmacfarlane.net/pandoc/
## or using your OS package manager (if running on Linux)

#install_required_packages()
#rsv.install.packages()

## loads dependency packages (which already must be installed)
load_required_packages()

library("BiocParallel")
register(SnowParam(threads))

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"g_test.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=rsv.vac.debug.mode)

evalsOptions("graph.output","svg")

set.seed(10)

rsv_run_project <- function() {
  
  report <<- PandocAT$new(author="tovchigrechkoa@medimmune.com",
                          title="RSV Surveilance Sequence Analysis",
                          incremental.save=F)
  
  
  rsv_process()
  report$save()
}
