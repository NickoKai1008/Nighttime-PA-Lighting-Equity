# V2 mean-model utilities. No legacy dispatcher or outcome fit runs on import.
suppressPackageStartupMessages({library(data.table); library(mgcv)})
V2_FIELDS <- c("road_mean","road_relative_iqr","walkability_global_z","NDVI",
  "ssx_choice_r5000","BD200","poi_density_500m","poi_type_mix_500m",
  "distance_to_branch_m","mapped_terminal_branch","isovist_fraction",
  "Pop_Female_20_50","gdp_2022","tavg","prcp","wspd","pres","rhum")
V2_WEATHER <- c("tavg","prcp","wspd","pres","rhum")
V2_DAY_ROOT <- "D:/codex\u8bca\u65ad/outputs/GAMM_day_night_spatial_20261028"
V2_LAND_ROOT <- "C:/Users/pc/.codex/visualizations/2026/10/23/01a167c1-3af6-7723-8e32-d4b8fbafb264/landuse_k4_k9"
v2_json <- function(path) jsonlite::fromJSON(path,simplifyVector=FALSE)
v2_sha <- function(path) {
  # Stream the actual bytes: NAS owner SID lookup can make file.access false
  # even when the authenticated connection can read the file normally.
  connection<-file(path,open="rb");on.exit(close(connection))
  unclass(as.character(openssl::sha256(connection)))
}
v2_record <- function(path) list(path=normalizePath(path,winslash="/",mustWork=TRUE),sha256=v2_sha(path))
v2_dir <- function(path) dir.create(path,recursive=TRUE,showWarnings=FALSE)
v2_retry <- function(fun) {
  for(delay in c(0,.15,.35,.75,1.5,3)) {
    if(delay>0)Sys.sleep(delay)
    if(isTRUE(suppressWarnings(fun())))return(TRUE)
  }
  FALSE
}
v2_commit <- function(tmp,path) {
  backup<-paste0(path,".previous-",Sys.getpid());had_old<-file.exists(path)
  if(had_old&&!v2_retry(function()file.rename(path,backup)))
    stop("Target is locked; old result preserved: ",path)
  if(!v2_retry(function()file.rename(tmp,path))) {
    if(had_old&&!file.exists(path))v2_retry(function()file.rename(backup,path))
    stop("Atomic result rename failed; temporary preserved: ",tmp)
  }
  if(had_old&&file.exists(backup))v2_retry(function()unlink(backup)==0L)
  invisible(path)
}
v2_json_write <- function(x,path) {
  v2_dir(dirname(path));tmp<-paste0(path,".partial-",Sys.getpid())
  jsonlite::write_json(x,tmp,auto_unbox=TRUE,pretty=TRUE,null="null",digits=17,na="null")
  v2_commit(tmp,path)
}
v2_rds_write <- function(x,path) {
  v2_dir(dirname(path));tmp<-paste0(path,".partial-",Sys.getpid());saveRDS(x,tmp,compress=FALSE)
  v2_commit(tmp,path)
}
v2_csv_write <- function(x,path) {
  v2_dir(dirname(path));tmp<-paste0(path,".partial-",Sys.getpid());fwrite(x,tmp,na="NA")
  v2_commit(tmp,path)
}
v2_id <- function(job) paste(job$scope,job$architecture,job$city,job$period,
  job$family,paste0("k",job$k),paste0("fold",job$fold),sprintf("r%04d",job$rep),sep="__")
v2_sources <- function(job) {
  base<-if(job$scope=="all25")V2_DAY_ROOT else file.path(V2_LAND_ROOT,job$scope)
  ip<-file.path(base,"inputs/INPUT_SPEC.json")
  mp<-file.path(base,job$period,"inputs/METADATA.rds")
  jp<-sub("\\.rds$",".json",mp)
  if(!all(file.exists(c(ip,mp,jp))))stop("Required V1 input/metadata absent.")
  spec<-v2_json(ip);mark<-v2_json(jp)
  if(!identical(spec$status,"complete")||
     !setequal(unlist(spec$active_covariates,use.names=FALSE),V2_FIELDS))stop("Input specification changed.")
  if(!identical(mark$status,"complete")||!identical(mark$metadata$sha256,v2_sha(mp)))
    stop("V1 metadata payload SHA mismatch.")
  # Verify every source identity claimed by the existing metadata marker.
  for(key in c("input","source_25_city_input","scaling_source")) {
    z<-mark[[key]]
    if(!is.null(z)&&!is.null(z$path)&&!identical(v2_sha(z$path),z$sha256))
      stop("V1 metadata source changed: ",key)
  }
  meta<-readRDS(mp)
  if(!all(V2_FIELDS%in%meta$scaling$field))stop("Required scaling columns absent.")
  declared<-sort(vapply(spec$cities,function(x)x$city,character(1)))
  # This execution ID is the same-cohort bridge. Expanded cohorts need a new variant ID.
  cohort<-sort(unlist(meta$cities,use.names=FALSE))
  if(!length(cohort)||!all(cohort%in%declared))stop("Invalid fixed bridge cohort.")
  if(job$architecture=="independent") {
    if(!job$city%in%declared)stop("Unknown independent city.")
    cohort<-intersect(cohort,job$city)
  }
  items<-Filter(function(x)x$city%in%cohort,spec$cities)
  list(spec=spec,items=items,declared=declared,cohort=cohort,
    scaling=as.data.table(meta$scaling)[field%in%V2_FIELDS],
    reference=list(path=mp,sha256=v2_sha(mp)),
    identity=list(input=v2_record(ip),metadata=v2_record(mp),metadata_marker=v2_record(jp)),
    supports=as.data.table(meta$supports))
}
v2_read_city <- function(item,period) {
  if(!identical(v2_sha(item$path),item$sha256))stop("Input SHA mismatch: ",item$city)
  keys<-c("city","grid_id","date","date_num","is_weekend","spatial_block_2km","spatial_block_5km")
  need<-unique(c(keys,setdiff(V2_FIELDS,V2_WEATHER),paste0(V2_WEATHER,"_",period),
    paste0(c("female","male"),"_sum_",period)))
  columns<-unlist(item$columns,use.names=FALSE)
  if(!length(columns)) {
    con<-gzfile(item$path,"rt");columns<-names(read.csv(con,nrows=0,check.names=FALSE));close(con)
  }
  if(!all(need%in%columns))stop("Missing input columns: ",paste(setdiff(need,columns),collapse=","))
  classes<-setNames(rep("NULL",length(columns)),columns)
  classes[need]<-"numeric";classes[intersect(need,c("city","grid_id","date","spatial_block_2km","spatial_block_5km"))]<-"character"
  con<-gzfile(item$path,"rt");on.exit(close(con))
  d<-as.data.table(read.csv(con,check.names=FALSE,stringsAsFactors=FALSE,colClasses=classes,
    na.strings=c("","NA","NaN")))
  if(nrow(d)!=item$rows||any(d$city!=item$city)||anyDuplicated(d[,.(grid_id,date)]))
    stop("Input count/key mismatch: ",item$city)
  for(f in V2_WEATHER)set(d,j=f,value=d[[paste0(f,"_",period)]])
  for(s in c("female","male"))setnames(d,paste0(s,"_sum_",period),paste0("y_",s))
  for(f in c(V2_FIELDS,"date_num","is_weekend","y_female","y_male"))
    if(any(!is.finite(d[[f]])))stop("Missing/nonfinite input; never NA-to-zero: ",item$city," ",f)
  if(any(d$y_female<0|d$y_male<0)||any(d$road_mean<0))stop("Negative PA/radiance.")
  for(b in c("spatial_block_2km","spatial_block_5km")) {
    if(anyNA(d[[b]])||any(!nzchar(d[[b]])))stop("Missing spatial block.")
    z<-d[[b]];prefix<-paste0(item$city,"::")
    if(any(!startsWith(z,prefix)))set(d,j=b,value=paste(item$city,z,sep="::"))
  }
  d[,c("date",paste0(V2_WEATHER,"_",period)):=NULL]
  d
}
v2_source_stat <- function(item) {
  f<-file.info(item$path)
  if(nrow(f)!=1L||is.na(f$size))stop("Input disappeared: ",item$path)
  list(city=item$city,path=normalizePath(item$path,winslash="/",mustWork=TRUE),
    expected_sha256=item$sha256,size=unname(f$size),mtime_numeric=unname(as.numeric(f$mtime)))
}
v2_raw_cached <- function(sources,job,build,log) {
  # Only immutable numeric/raw columns produced by v2_read_city: no z scaling,
  # fold/outcome filtering, FRW weights, or fitted quantities enter this cache.
  identity<-list(input=sources$identity$input,period=job$period,fields=as.list(V2_FIELDS),
    cohort=as.list(sources$cohort),payloads=lapply(sources$items,function(x)
      list(city=x$city,path=x$path,sha256=x$sha256,rows=x$rows)),
    read_city_body_sha256=digest::digest(body(v2_read_city),algo="sha256"))
  fingerprint<-digest::digest(identity,algo="sha256")
  cache_label<-paste(job$scope,job$period,job$city,sep="__")
  if(!is.null(job$variant)&&job$variant=="expanded")cache_label<-paste0(cache_label,"__variant_expanded")
  cache_dir<-file.path(job$root,"cache","raw",cache_label)
  path<-file.path(cache_dir,"RAW_WIDE.rds");marker<-file.path(cache_dir,"RAW_WIDE.json")
  current_stats<-lapply(sources$items,v2_source_stat)
  if(file.exists(path)&&file.exists(marker)) {
    prior<-v2_json(marker)
    # Reusing a prior byte-verified input requires unchanged expected hashes AND
    # size/mtime. A change forces full source hashing through the ordinary reader.
    stat_ok<-identical(digest::digest(current_stats,algo="sha256"),prior$source_stat_fingerprint)
    if(identical(prior$status,"complete")&&identical(prior$fingerprint,fingerprint)&&stat_ok&&
       identical(v2_sha(path),prior$cache$sha256)) {
      log("RAW_CACHE_HIT",path);d<-readRDS(path);setDT(d)
      if(nrow(d)!=sum(vapply(sources$items,function(x)as.numeric(x$rows),numeric(1)))||
         any(startsWith(names(d),"z_"))||any(c("fold_id","case_weight")%in%names(d)))
        stop("Raw cache schema/count mismatch.")
      return(list(data=d,identity=list(marker=v2_record(marker),cache=prior$cache,
        input_payload_verification="Previously SHA-verified bytes; expected SHA, size and mtime unchanged; cached RDS SHA checked now.")))
    }
    log("RAW_CACHE_STALE_REBUILD",path)
  }
  d<-build(sources)
  if(any(startsWith(names(d),"z_"))||any(c("fold_id","case_weight")%in%names(d)))stop("Raw cache contains derived model columns.")
  after_stats<-lapply(sources$items,v2_source_stat)
  if(!identical(current_stats,after_stats))stop("Source changed while raw cache was built.")
  v2_rds_write(d,path)
  mark<-list(status="complete",fingerprint=fingerprint,identity=identity,
    source_stat_fingerprint=digest::digest(after_stats,algo="sha256"),source_stats=after_stats,
    cache=v2_record(path),rows=nrow(d),columns=as.list(names(d)),
    source_payloads_byte_verified_by="v2_read_city sha256 checks before parsing each source CSV")
  v2_json_write(mark,marker);log("RAW_CACHE_CREATED",path)
  list(data=d,identity=list(marker=v2_record(marker),cache=mark$cache,
    input_payload_verification="All source payload SHA-256 values checked on raw-cache creation."))
}
v2_attach_folds <- function(d,root,required) {
  path<-file.path(root,"inputs/FOLD_MAP.csv")
  if(!file.exists(path)) {
    if(required)stop("FOLD_MAP.csv required for CV/paired FRW.")
    d[,fold_id:=NA_integer_];return(list(data=d,identity=NULL,blocks=sort(unique(d$spatial_block_2km))))
  }
  f<-fread(path,colClasses=c(city="character",grid_id="character",spatial_block_2km="character",spatial_block_5km="character"))
  if(!all(c("city","grid_id","fold_id","spatial_block_2km")%in%names(f))||anyDuplicated(f[,.(city,grid_id)]))
    stop("Invalid fold map keys.")
  f[,fold_id:=as.integer(fold_id)]
  if(any(!is.na(f$fold_id)&!f$fold_id%in%1:5))stop("Fold ID outside 1..5.")
  token<-function(value,city) {
    # Legacy tables can carry city:i:j, city::i:j, or city::city:i:j.
    # Remove only repeated copies of the row's own city, preserving spatial indices.
    value<-as.character(value)
    for(pass in 1:8) {
      hit<-startsWith(value,paste0(city,":"))
      if(!any(hit))break
      value[hit]<-sub("^:+","",substring(value[hit],nchar(city[hit])+1L))
    }
    if(anyNA(value)||any(!nzchar(value)))stop("Missing/empty block spatial token.")
    value
  }
  map_row<-f[d,on=.(city,grid_id),which=TRUE]
  if(anyNA(map_row))stop("Observed grid absent from the canonical fold/block registry.")
  for(b in c("spatial_block_2km","spatial_block_5km")) {
    if(!b%in%names(f))stop("Fold registry missing canonical block: ",b)
    original_token<-token(d[[b]],d$city);canonical_token<-token(f[[b]],f$city)
    if(any(original_token!=canonical_token[map_row]))
      stop("Source/fold spatial indices disagree for ",b,"; do not repair by label alone.")
    # One registry row per city-grid plus equality on EVERY date verifies that all
    # repeated grid-days keep identical physical block indices before canonicalizing.
    set(f,j=b,value=paste(f$city,canonical_token,sep="::"))
    set(d,j=b,value=f[[b]][map_row])
  }
  set(d,j="fold_id",value=f$fold_id[map_row])
  if(required&&anyNA(d$fold_id))stop("Missing fold assignment in requested training/validation cohort.")
  # Pair all scopes, sexes and dates by one globally ordered block registry.
  list(data=d,identity=v2_record(path),blocks=sort(unique(f$spatial_block_2km)))
}
v2_scale <- function(d,scaling) {
  for(j in seq_len(nrow(scaling))) {
    r<-scaling[j]
    if(!is.finite(r$sample_sd)||r$sample_sd<=0)stop("Invalid scale: ",r$field)
    set(d,j=paste0("z_",r$field),value=(d[[r$field]]-r$center)/r$sample_sd)
  }
  d
}
v2_training_scaling <- function(d,allow_constant=FALSE) rbindlist(lapply(V2_FIELDS,function(f) {
  x<-d[[f]];s<-sd(x);constant<-is.finite(s)&&s==0
  # A constant independent-city control is registered and removed before fitting.
  if(constant&&allow_constant)s<-1
  data.table(field=f,center=mean(x),sample_sd=s,N_pairs=length(x),constant_training_field=constant)
}))
v2_weights <- function(d,rep,blocks) {
  if(rep==0L)return(rep.int(1,nrow(d)))
  set.seed(20261102L+as.integer(rep));w<-rgamma(length(blocks),shape=1,rate=1);w<-w/mean(w)
  index<-match(d$spatial_block_2km,blocks)
  if(anyNA(index))stop("FRW block missing from the global registry.")
  w[index]
}
v2_formula <- function(architecture,k,controls,ncity) {
  baseline<-if(architecture=="city_deviation"||ncity==1L)"1"else"0+city_fac"
  terms<-c(baseline,controls,sprintf("s(z_road_mean,k=%d,bs='cr')",k),"s(date_num,k=5,bs='cr')")
  if(architecture=="city_deviation") {
    if(ncity<2L)stop("City deviation requires at least two estimable cities.")
    terms<-c(terms,sprintf("s(city_fac,z_road_mean,k=%d,bs='sz',xt=list(bs='cr'),id='dev')",k))
  }
  as.formula(paste("y~",paste(terms,collapse="+")),env=asNamespace("mgcv"))
}
v2_rank <- function(d,controls,architecture,ncity) {
  terms<-c(if(architecture=="city_deviation"||ncity==1L)"1"else"0+city_fac",controls)
  f<-as.formula(paste("~",paste(terms,collapse="+")));R<-NULL
  for(first in seq.int(1L,nrow(d),by=5000L)) {
    X<-model.matrix(f,data=as.data.frame(d[seq.int(first,min(nrow(d),first+4999L))]))
    qrX<-qr(if(is.null(R))X else rbind(R,X),LAPACK=TRUE)
    R<-qr.R(qrX)[,order(qrX$pivot),drop=FALSE]
  }
  z<-qr(R,tol=1e-9,LAPACK=FALSE);p<-ncol(R)
  list(rank=z$rank,p=p,aliases=if(z$rank<p)colnames(R)[z$pivot[seq.int(z$rank+1L,p)]]else character())
}
v2_poisson_deviance <- function(y,mu) {
  if(any(!is.finite(mu)|mu<0))stop("Invalid mean prediction.")
  out<-2*(mu-y);positive<-y>0
  out[positive]<-2*(y[positive]*log(y[positive]/mu[positive])-y[positive]+mu[positive])
  pmax(out,0)
}
v2_compact <- function(object,probe) {
  before<-predict(object,probe,type="lpmatrix",discrete=FALSE)
  object$model<-object$model[0,,drop=FALSE]
  for(f in c("fitted.values","linear.predictors","residuals","prior.weights","weights","y",
             "AR.start","dinfo","X","hat","offset"))object[[f]]<-NULL
  after<-predict(object,probe,type="lpmatrix",discrete=FALSE)
  if(!isTRUE(all.equal(before,after,tolerance=1e-12)))stop("Compact model changed prediction design.")
  object
}
v2_grid_residuals <- function(object,keys,sex,period) {
  y<-object$y;mu<-object$fitted.values
  phi<-as.numeric(object$sig2)
  if(length(phi)!=1L||!is.finite(phi)||phi<=0)stop("Invalid fitted dispersion for Pearson residuals.")
  variance<-object$family$variance(mu)
  if(length(variance)!=length(y)||any(!is.finite(variance)|variance<0))stop("Invalid fitted family variance.")
  pearson<-(y-mu)/sqrt(phi*variance)
  pearson[!is.finite(pearson)]<-NA_real_
  z<-data.table(city=as.character(keys$city),grid_id=as.character(keys$grid_id),
    date_num=keys$date_num,raw_resid=y-mu,pearson=pearson,y=y,mu=mu)
  ans<-z[,.(n_rows=.N,n_dates=uniqueN(date_num),mean_raw_resid=mean(raw_resid),
    mean_Pearson_resid=if(all(is.na(pearson)))NA_real_ else mean(pearson,na.rm=TRUE),
    n_Pearson_finite=sum(is.finite(pearson)),observed_mean=mean(y),predicted_mean=mean(mu)),
    by=.(city,grid_id)]
  ans[,`:=`(sex=sex,period=period,Pearson_dispersion=phi,
    residual_status=ifelse(n_Pearson_finite==n_rows,"fitted_grid_mean","some_Pearson_residuals_undefined"))]
  ans
}
v2_restore_sample <- function(object,data,indices,formula,object_indices=indices,source_n=nrow(data)) {
  z<-object
  gp<-get("interpret.gam",asNamespace("mgcv"))(formula)
  mf<-model.frame(gp$fake.formula,data=as.data.frame(data[indices]),
    na.action=na.fail,drop.unused.levels=FALSE)
  attr(mf,"terms")<-object$terms;attr(mf,"na.action")<-NULL;z$model<-mf
  for(f in c("y","fitted.values","linear.predictors","prior.weights","weights","residuals","offset"))
    if(!is.null(object[[f]])&&length(object[[f]])==source_n)z[[f]]<-object[[f]][object_indices]
  z$na.action<-NULL;z
}
v2_subsample_diagnostics <- function(object,data,formula,sex,period,prepared_indices=NULL,full_training_rows=nrow(data)) {
  n<-as.integer(full_training_rows)
  if(is.null(prepared_indices)) {
    set.seed(20261102L+match(sex,c("female","male"))+10L*match(period,c("day","night")))
    indices<-if(n<=5000L)seq_len(n)else sort(sample.int(n,5000L))
    object_indices<-indices
  }else {
    if(length(prepared_indices)!=nrow(data))stop("Prepared diagnostic sample/index mismatch.")
    indices<-seq_len(nrow(data));object_indices<-as.integer(prepared_indices)
  }
  z<-v2_restore_sample(object,data,indices,formula,object_indices=object_indices,source_n=n)
  errors<-list();kc<-tryCatch(mgcv::k.check(z,subsample=5000L,n.rep=100L),error=function(e) {
    errors$k_check<<-conditionMessage(e);NULL
  })
  cv<-tryCatch(mgcv::concurvity(z,full=TRUE),error=function(e) {
    errors$concurvity<<-conditionMessage(e);NULL
  })
  ktable<-if(is.null(kc))data.table()else data.table(sex=sex,smooth=rownames(kc),
    k_prime=kc[,1L],edf=kc[,2L],k_index=kc[,3L],p_value=kc[,4L],
    diagnostic_rows=length(indices),full_training_rows=n,n_rep=100L,
    method="Fixed-seed at-most-5000 training-row residual diagnostic; factor smooth k checks may be NA")
  ctable<-if(is.null(cv))data.table()else data.table(sex=sex,term=colnames(cv),
    worst=cv["worst",],observed=cv["observed",],estimate=cv["estimate",],
    diagnostic_rows=length(indices),full_training_rows=n,
    method="Concurvity on fixed-seed at-most-5000-row training model frame; association overlap, not importance")
  list(status=if(length(errors))"partial_or_failed"else"complete_subsample_diagnostic",
    k_check=ktable,concurvity=ctable,errors=errors,
    sample_rows=length(indices),full_rows=n,
    note="These model diagnostics do not treat grid-date rows as independent spatial sites and are not formal spatial tests.")
}
v2_smoke_equivalence <- function(object,compact,data,formula,knots,job,threads) {
  # Run only in --op smoke on synthetic data, never repeat the full production fit.
  fam<-if(job$family=="quasipoisson_log")quasipoisson(link="log")else mgcv::tw(link="log")
  direct<-mgcv::bam(formula,data=data,weights=case_weight,family=fam,method="fREML",
    discrete=10000L,nthreads=threads,gc.level=2L,chunk.size=5000L,
    na.action=na.fail,drop.unused.levels=FALSE,knots=list(z_road_mean=knots))
  index<-unique(as.integer(round(seq(1,nrow(data),length.out=min(64L,nrow(data))))))
  probe<-data[index]
  dp<-predict(direct,probe,type="response",discrete=FALSE)
  sp<-predict(object,probe,type="response",discrete=FALSE)
  cp<-predict(compact,probe,type="response",discrete=FALSE)
  dl<-predict(direct,probe,type="link",discrete=FALSE)
  sl<-predict(object,probe,type="link",discrete=FALSE)
  cd<-predict(compact,probe,type="lpmatrix",discrete=FALSE)
  sd<-predict(object,probe,type="lpmatrix",discrete=FALSE)
  response_error<-max(abs(dp-sp));compact_error<-max(abs(sp-cp));link_error<-max(abs(dl-sl))
  design_error<-max(abs(cd-sd));tol<-1e-7
  dev_ok<-isTRUE(all.equal(as.numeric(direct$deviance),as.numeric(object$deviance),tolerance=tol))
  null_ok<-isTRUE(all.equal(as.numeric(direct$null.deviance),as.numeric(object$null.deviance),tolerance=tol))
  fam_ok<-identical(direct$family$family,object$family$family)&&identical(direct$family$link,object$family$link)
  ok<-isTRUE(direct$converged)&&isTRUE(object$converged)&&fam_ok&&dev_ok&&null_ok&&
    response_error<=tol*max(1,max(abs(dp)))&&link_error<=tol*max(1,max(abs(dl)))&&
    compact_error<=1e-12*max(1,max(abs(sp)))&&design_error<=1e-12
  list(status=if(ok)"passed"else"failed",synthetic_data=TRUE,N=nrow(data),probe_rows=length(index),
    family=object$family$family,architecture=job$architecture,k=job$k,
    max_abs_full_vs_split_response=response_error,max_abs_full_vs_split_link=link_error,
    max_abs_split_vs_compact_response=compact_error,max_abs_split_vs_compact_lpmatrix=design_error,
    direct_deviance=direct$deviance,split_deviance=object$deviance,
    direct_null_deviance=direct$null.deviance,split_null_deviance=object$null.deviance,
    family_equal=fam_ok,deviance_equal=dev_ok,null_deviance_equal=null_ok,
    full_vs_split_relative_tolerance=tol,compact_tolerance=1e-12,
    comparison="Direct full-frame bam versus fit=FALSE G + response-only mf fit + compact prediction")
}
