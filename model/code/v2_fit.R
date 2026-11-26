# Direct nonnegative mean GAMs. Invoke under the shared Windows Job/memory guard.
argv<-commandArgs(trailingOnly=TRUE)
script_arg<-grep("^--file=",commandArgs(),value=TRUE)
if(!length(script_arg))stop("Invoke with Rscript.")
SCRIPT<-normalizePath(sub("^--file=","",script_arg[[1L]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(SCRIPT),"v2_model_common.R"),encoding="UTF-8")
source(file.path(dirname(SCRIPT),"v2_expanded_sources.R"),encoding="UTF-8")
options(stringsAsFactors=FALSE)
threads<-as.integer(Sys.getenv("V2_NTHREADS","2"))
if(!is.finite(threads)||!threads%in%1:4)stop("V2_NTHREADS must be 1..4.")
setDTthreads(1L)
Sys.setenv(OMP_NUM_THREADS=as.character(threads),OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
V2_SOLVER_REPAIR<-list(version="strict_bgam_convergence_and_intercept_retry_v1",
  predecessor_engine_sha256="b1e1845c781f1b8bdaa70a7a0e71a3e5b248a620344964d15032372f9ff7d2d5",
  unchanged="Input rows, encoded-zero response, controls, formula, knots, scaling, family, weights and fREML objective",
  initial_attempt="Unchanged default initialization, epsilon=1e-7, maxit=200",
  retry="Only after solver failure: weighted city-mean intercept starts (global mean when no city FE), other coefficients zero; maxit=400 then 1000; epsilon remains 1e-7; fresh family, all sp/p re-estimated",
  acceptance="converged AND mgcv.conv AND finite coefficient/mean/Vp/deviance AND no solver-failure warning")

parse_job <- function(args) {
  if(length(args)%%2L)stop("Every CLI option needs a value.")
  keys<-sub("^--","",args[seq.int(1L,length(args),by=2L)])
  if(anyDuplicated(keys)||any(!keys%in%c("op","scope","architecture","city","period","family","k","fold","rep","root","variant")))
    stop("Unknown/repeated CLI option.")
  a<-setNames(as.list(args[seq.int(2L,length(args),by=2L)]),keys)
  required<-c("op","scope","architecture","city","period","family","k","fold","rep","root")
  if(!all(required%in%names(a)))stop("Missing CLI option: ",paste(setdiff(required,names(a)),collapse=","))
  a$variant<-v2_job_variant(a)
  a$k<-as.integer(a$k);a$rep<-as.integer(a$rep)
  if(!a$op%in%c("smoke","fit")||!a$scope%in%c("all25","commercial","residential","recreational")||
     !a$architecture%in%c("shared","city_deviation","independent")||!a$period%in%c("day","night")||
     !a$family%in%c("quasipoisson_log","tweedie_log")||!a$k%in%4:9||
     !a$fold%in%c("all",as.character(1:5))||!is.finite(a$rep)||a$rep<0L)
    stop("Invalid CLI value.")
  if((a$architecture=="independent")==(a$city=="pooled"))stop("Independent requires a city; other architectures require pooled.")
  a$root<-normalizePath(a$root,winslash="/",mustWork=TRUE);a
}
job<-parse_job(argv);ID<-v2_variant_id(job)
OUT<-file.path(job$root,if(job$op=="smoke")"smoke"else"models",ID);v2_dir(OUT)
LOG<-file.path(OUT,"RUN.log")
log_stage<-function(...) {
  line<-paste(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),paste(...,collapse=" "))
  cat(line,"\n");cat(line,"\n",file=LOG,append=TRUE);flush.console()
}
status_write<-function(status,details=list())v2_json_write(c(list(status=status,task_id=ID,job=job,
  updated_local=format(Sys.time(),"%Y-%m-%d %H:%M:%S")),details),file.path(OUT,"STATUS.json"))

smoke_data <- function(k) {
  set.seed(20261102L);n<-2400L
  d<-data.table(city=rep(c("smoke_a","smoke_b","smoke_c"),each=n/3),grid_id=as.character(seq_len(n)),
    date_num=rep(1:20,length.out=n),is_weekend=rep(c(0,0,0,0,0,1,1),length.out=n))
  for(f in V2_FIELDS)set(d,j=f,value=rnorm(n))
  d[,road_mean:=exp(road_mean-7)];d[,road_relative_iqr:=abs(road_relative_iqr)]
  d[,spatial_block_2km:=paste(city,rep(1:40,length.out=.N),sep="::")]
  d[,spatial_block_5km:=paste(city,rep(1:20,length.out=.N),sep="::")]
  d[,`:=`(y_female=rpois(.N,exp(-1.6+0.3*log(road_mean/0.001))),
          y_male=rpois(.N,exp(0.2+0.2*log(road_mean/0.001))),fold_id=rep(1:5,length.out=.N))]
  # Fractional response tests the archived-length allocation estimand explicitly.
  d[,`:=`(y_female=y_female*0.37,y_male=y_male*0.37)]
  if(job$architecture=="independent")d<-d[city=="smoke_a"]
  d
}
make_data <- function(sources) {
  parts<-vector("list",length(sources$items))
  for(i in seq_along(sources$items)) {
    item<-sources$items[[i]];log_stage("READ",item$city,item$rows)
    parts[[i]]<-v2_read_city(item,job$period)
  }
  d<-rbindlist(parts,use.names=TRUE);rm(parts);gc();d
}
prepare_sex <- function(wide,sex) {
  controls<-c(paste0("z_",setdiff(V2_FIELDS,"road_mean")),"is_weekend")
  cols<-c("city","grid_id","spatial_block_2km","date_num","fold_id",controls,"z_road_mean")
  d<-copy(wide[,cols,with=FALSE]);d[,y:=wide[[paste0("y_",sex)]]]
  d[,case_weight:=wide$case_weight]
  train<-if(job$fold=="all")rep.int(TRUE,nrow(d))else d$fold_id!=as.integer(job$fold)
  tab<-d[,.(n=.N,n_positive=sum(y>0),PA_sum=sum(y),unique_L=uniqueN(z_road_mean),
    blocks=uniqueN(spatial_block_2km)),by=.(city,training=train)]
  trainstats<-tab[training==TRUE]
  trainstats[,status:=ifelse(n_positive==0,"all_zero_FE_boundary","estimable_mean_cell")]
  tab[,status:=ifelse(training,ifelse(n_positive==0,"all_zero_FE_boundary","estimable_mean_cell"),"validation_rows")]
  keepcity<-trainstats$city[trainstats$n_positive>0]
  d[,city_fac:=factor(city,levels=sort(keepcity))]
  # Keep held-out rows in the wide object for honest boundary/unseen-city scores.
  tr<-d[train&city%in%keepcity]
  dropped<-data.table(field=character(),reason=character())
  if(!nrow(tr))return(list(data=tr,training=tab,status="no_estimable_city_sex",controls=controls,dropped=dropped))
  if(job$architecture=="independent") {
    constants<-controls[vapply(controls,function(f)uniqueN(tr[[f]])<2L,logical(1))]
    if(length(constants))dropped<-rbind(dropped,data.table(field=constants,reason="constant_in_training_city"))
    controls<-setdiff(controls,constants)
  }
  if(uniqueN(tr$z_road_mean)<job$k||uniqueN(tr$date_num)<5L)
    return(list(data=tr,training=tab,status="insufficient_L_or_date_basis_support",controls=controls,dropped=dropped))
  rank<-v2_rank(tr,controls,job$architecture,length(keepcity))
  if(rank$rank<rank$p&&job$architecture=="independent") {
    aliases<-intersect(rank$aliases,controls)
    if(length(aliases)) {
      dropped<-rbind(dropped,data.table(field=aliases,reason="exact_parametric_linear_dependence_fixed_order"))
      controls<-setdiff(controls,aliases);rank<-v2_rank(tr,controls,job$architecture,length(keepcity))
    }
  }
  list(data=tr,training=tab,status=if(rank$rank==rank$p)"ready"else"parametric_rank_deficient",
    controls=controls,dropped=dropped,rank=rank,city_levels=sort(keepcity))
}
fit_sex <- function(prepared,sex,knots) {
  d<-prepared$data
  model_n<-nrow(d)
  form<-v2_formula(job$architecture,job$k,prepared$controls,length(prepared$city_levels))
  fam<-if(job$family=="quasipoisson_log")quasipoisson(link="log")else mgcv::tw(link="log")
  warn<-character();started<-Sys.time();log_stage("SETUP",sex,nrow(d),deparse(form,width.cutoff=500L))
  G<-withCallingHandlers(mgcv::bam(form,data=d,weights=case_weight,family=fam,method="fREML",
    discrete=10000L,nthreads=threads,gc.level=2L,chunk.size=5000L,na.action=na.fail,
    drop.unused.levels=FALSE,knots=list(z_road_mean=knots),fit=FALSE),
    warning=function(w){warn<<-unique(c(warn,conditionMessage(w)));invokeRestart("muffleWarning")})
  probe<-d[seq_len(min(8L,nrow(d)))]
  city_codes<-as.integer(d$city_fac);city_levels<-levels(d$city_fac)
  requested_weights<-d$case_weight
  residual_keys<-diagnostic_data<-diagnostic_indices<-NULL
  if(job$fold=="all"&&job$rep==0L) {
    residual_keys<-d[,.(city,grid_id,date_num)]
    set.seed(20261102L+match(sex,c("female","male"))+10L*match(job$period,c("day","night")))
    diagnostic_indices<-if(model_n<=5000L)seq_len(model_n)else sort(sample.int(model_n,5000L))
    diagnostic_data<-d[diagnostic_indices]
  }
  # rho=0 fit uses response-only model frame; the compressed design remains in G.
  if(!isTRUE(G$discretize)||!isTRUE(all.equal(as.numeric(G$mf$y),as.numeric(G$y),tolerance=0)))
    stop("Unexpected discrete setup response.")
  slim<-G$mf["y"];attr(slim,"terms")<-attr(G$mf,"terms");attr(slim,"na.action")<-attr(G$mf,"na.action")
  G$mf<-slim
  # prepared is an environment owned by the caller. Remove the actual full frame,
  # not just this function's local alias. Only synthetic smoke retains it for the
  # direct full-frame equivalence comparison; production keeps tiny diagnostics.
  if(job$op!="smoke")prepared$data<-NULL
  rm(slim,d);gc();log_stage("FIT",sex)
  setup_warnings<-warn;attempts<-list();object<-NULL;family_template<-serialize(G$family,NULL)
  solver_warning<-function(x)any(grepl("algorithm did not converge|not converg|step failure|iteration limit|non-finite coefficients|outer.*fail",
    x,ignore.case=TRUE))
  for(attempt in 1:3) {
    attempt_warnings<-character();attempt_error<-NULL;maximum<-c(200L,400L,1000L)[[attempt]]
    if(attempt==1L)initial_coef<-NULL else {
      # This changes the start only. Zero-boundary city-sex cells have already
      # been registered and excluded, so each estimated intercept has positive PA.
      initial_coef<-numeric(ncol(G$X));cn<-colnames(G$X)
      numerator<-rowsum(as.numeric(G$y)*as.numeric(G$w),city_codes,reorder=FALSE)
      denominator<-rowsum(as.numeric(G$w),city_codes,reorder=FALSE)
      means<-as.numeric(numerator/denominator)
      if(any(!is.finite(means)|means<=0))stop("Invalid weighted intercept initialization: ",sex)
      if(job$architecture!="city_deviation"&&length(city_levels)>1L) {
        columns<-match(paste0("city_fac",city_levels[as.integer(rownames(numerator))]),cn)
        if(anyNA(columns))stop("City intercept initialization columns absent.")
        initial_coef[columns]<-log(means)
      }else {
        column<-match("(Intercept)",cn)
        if(is.na(column))stop("Global intercept initialization column absent.")
        initial_coef[[column]]<-log(sum(as.numeric(G$y)*as.numeric(G$w))/sum(G$w))
      }
      rm(numerator,denominator,means);gc()
    }
    # Keep family functions augmented during bam setup, but copy the pending
    # theta environment afresh so no failed estimate enters the next attempt.
    G$family<-unserialize(family_template)
    if(any(G$sp>=0)||!identical(as.numeric(G$w),as.numeric(requested_weights)))
      stop("Retry setup no longer has pending sp and the original prior weights.")
    log_stage("SOLVER_ATTEMPT",sex,attempt,"maxit",maximum,
      if(attempt==1L)"default_start"else"weighted_intercept_start")
    attempt_started<-Sys.time()
    object<-tryCatch(withCallingHandlers(mgcv::bam(G=G,family=G$family,coef=initial_coef,
      method="fREML",nthreads=threads,gc.level=2L,chunk.size=5000L,rho=0,
      control=mgcv::gam.control(trace=TRUE,epsilon=1e-7,maxit=maximum)),
      warning=function(w){attempt_warnings<<-unique(c(attempt_warnings,conditionMessage(w)));invokeRestart("muffleWarning")}),
      error=function(e){attempt_error<<-conditionMessage(e);NULL})
    finite_fit<-!is.null(object)&&all(is.finite(coef(object)))&&all(is.finite(object$Vp))&&
      all(is.finite(object$fitted.values)&object$fitted.values>=0)&&is.finite(object$deviance)&&
      all(is.finite(object$sp)&object$sp>=0)
    accepted<-finite_fit&&isTRUE(object$converged)&&isTRUE(object$mgcv.conv)&&!solver_warning(attempt_warnings)
    attempts[[attempt]]<-list(attempt=attempt,start=if(attempt==1L)"default"else"weighted_intercepts_only",
      maxit=maximum,epsilon=1e-7,accepted=accepted,error=attempt_error,warnings=as.list(attempt_warnings),
      converged=if(!is.null(object))isTRUE(object$converged)else FALSE,
      mgcv_conv=if(!is.null(object))isTRUE(object$mgcv.conv)else FALSE,
      iterations=if(!is.null(object))object$iter else NULL,finite_fit=finite_fit,
      maximum_outer_gradient=if(!is.null(object)&&length(object$outer.info$grad))max(abs(object$outer.info$grad))else NULL,
      elapsed_seconds=as.numeric(difftime(Sys.time(),attempt_started,units="secs")))
    if(accepted){warn<-unique(c(setup_warnings,attempt_warnings));break}
    if(!is.null(attempt_error)&&grepl("alloc|R_Calloc|memory",attempt_error,ignore.case=TRUE))
      stop("Memory allocation failed; no repeated optimizer launch: ",attempt_error)
    rm(object);object<-NULL;gc()
  }
  rm(G,initial_coef,family_template);gc()
  if(is.null(object)) {
    v2_json_write(list(status="failed_solver",sex=sex,solver_repair=V2_SOLVER_REPAIR,attempts=attempts),
      file.path(OUT,paste0("SOLVER_FAILURE_",toupper(sex),".json")))
    stop("Mean fit failed strict convergence checks after controlled retries: ",sex)
  }
  if(!identical(as.numeric(object$prior.weights),as.numeric(requested_weights)))
    stop("Fitted prior weights differ from requested paired row weights: ",sex)
  # Explicit family verification avoids the legacy bam(G=G) family-summary pitfall.
  if(job$family=="quasipoisson_log"&&!identical(object$family$family,"quasipoisson"))stop("Wrong quasi family.")
  if(job$family=="tweedie_log"&&!grepl("Tweedie",object$family$family,fixed=TRUE))stop("Wrong Tweedie family.")
  y<-object$y;mu<-object$fitted.values
  metric<-data.table(city=city_levels[city_codes],y=y,mu=mu,dev=v2_poisson_deviance(y,mu))[
    ,.(n=.N,n_positive=sum(y>0),observed_sum=sum(y),predicted_sum=sum(mu),
       poisson_deviance_sum=sum(dev),poisson_deviance_mean=mean(dev),
       squared_error_sum=sum((y-mu)^2),absolute_error_sum=sum(abs(y-mu))),by=city]
  metric[,`:=`(sex=sex,period=job$period,dataset="training",task_id=ID)]
  edf<-data.table(term=names(object$edf),estimate=unname(coef(object)),edf=unname(object$edf),
    model_based_se=sqrt(pmax(diag(object$Vp),0)),sex=sex)
  summary<-list(status="fitted",N=length(y),positive=sum(y>0),encoded_zero=sum(y==0),
    family=object$family$family,link=object$family$link,deviance=object$deviance,
    null_deviance=object$null.deviance,explained_deviance=1-object$deviance/object$null.deviance,
    EDF=sum(object$edf),scale=object$sig2,sp=as.list(object$sp),warnings=as.list(warn),
    p=if(job$family=="tweedie_log")object$family$getTheta(TRUE)else NULL,
    elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),
    formula=paste(deparse(form,width.cutoff=500L),collapse=" "),
    controls=as.list(prepared$controls),city_levels=as.list(prepared$city_levels),rank=prepared$rank,
    prior_weights_verified=TRUE,prior_weight_sum=sum(object$prior.weights),
    prior_weight_sha256=digest::digest(as.numeric(object$prior.weights),algo="sha256"),
    mgcv_conv=isTRUE(object$mgcv.conv),iterations=object$iter,
    solver_repair=V2_SOLVER_REPAIR,solver_attempts=attempts,solver_recovered=length(attempts)>1L,
    solver_status=list(status="strict_checks_passed",converged=isTRUE(object$converged),object_converged=isTRUE(object$converged),
      mgcv_conv=isTRUE(object$mgcv.conv),finite_coefficient_mean_covariance_deviance=TRUE,
      final_attempt_solver_warning=FALSE,attempt_count=length(attempts),
      objective_unchanged=TRUE,epsilon=1e-7,sp_reestimated=TRUE,Tweedie_p_reestimated=job$family=="tweedie_log",
      coefficients=list(available=TRUE,count=length(coef(object)),nonfinite=sum(!is.finite(coef(object)))),
      Vp=list(available=TRUE,rows=nrow(object$Vp),columns=ncol(object$Vp),nonfinite=sum(!is.finite(object$Vp))),
      negative_Vp_diagonal=sum(diag(object$Vp)<0),minimum_Vp_diagonal=min(diag(object$Vp)),
      smoothing_parameters=list(available=TRUE,count=length(object$sp),nonfinite=sum(!is.finite(object$sp)),negative=sum(object$sp<0)),
      outer_info=list(available=!is.null(object$outer.info),
        maximum_absolute_gradient=if(length(object$outer.info$grad))max(abs(object$outer.info$grad))else NULL,
        nonfinite_gradient=if(length(object$outer.info$grad))sum(!is.finite(object$outer.info$grad))else NULL)),
    uncertainty="Model-based coefficient covariance only; paired complete-refit FRW is a separate stage.")
  grid_residuals<-data.table();diagnostics<-NULL
  if(job$fold=="all"&&job$rep==0L) {
    log_stage("GRID_RESIDUALS",sex)
    grid_error<-NULL
    grid_residuals<-tryCatch(v2_grid_residuals(object,residual_keys,sex,job$period),error=function(e) {
      grid_error<<-conditionMessage(e);data.table()
    })
    log_stage("SUBSAMPLE_DIAGNOSTICS",sex)
    diagnostics<-tryCatch(v2_subsample_diagnostics(object,diagnostic_data,form,sex,job$period,
      prepared_indices=diagnostic_indices,full_training_rows=model_n),error=function(e)
      list(status="failed",errors=list(sample_restore=conditionMessage(e)),k_check=data.table(),concurvity=data.table()))
    summary$grid_residuals<-list(status=if(is.null(grid_error))"complete"else"failed",error=grid_error,grids=nrow(grid_residuals),
      undefined_Pearson_rows=if(nrow(grid_residuals))sum(grid_residuals$n_rows-grid_residuals$n_Pearson_finite)else NULL,
      method="Residuals aggregated to city-grid means before spatial diagnostics; fitted family variance times estimated dispersion.")
    summary$subsample_diagnostics<-diagnostics[setdiff(names(diagnostics),c("k_check","concurvity"))]
  }else {
    summary$grid_residuals<-list(status="not_requested_for_CV_or_bootstrap")
    summary$subsample_diagnostics<-list(status="not_requested_for_CV_or_bootstrap")
  }
  # Reconstruct only a tiny valid model frame for post-fit lpmatrix prediction.
  gp<-get("interpret.gam",asNamespace("mgcv"))(form)
  mf<-model.frame(gp$fake.formula,data=as.data.frame(probe),na.action=na.fail,drop.unused.levels=FALSE)
  attr(mf,"terms")<-object$terms;object$model<-mf
  compact<-v2_compact(object,probe)
  if(job$op=="smoke") {
    log_stage("NUMERICAL_EQUIVALENCE",sex)
    qa<-v2_smoke_equivalence(object,compact,prepared$data,form,knots,job,threads)
    v2_json_write(qa,file.path(OUT,paste0("NUMERICAL_QA_",toupper(sex),".json")))
    summary$numerical_qa<-qa
    if(!identical(qa$status,"passed"))stop("Full/split/compact numerical equivalence failed: ",sex)
  }else summary$numerical_qa<-list(status="validated_in_separate_synthetic_smoke_not_repeated_in_production")
  log_stage("FIT_DONE",sex,summary$elapsed_seconds)
  list(model=compact,summary=summary,coef=edf,training_score=metric,
    grid_residuals=grid_residuals,diagnostics=diagnostics,solver_attempts=rbindlist(lapply(attempts,function(a)
      data.table(sex=sex,attempt=a$attempt,start=a$start,maxit=a$maxit,epsilon=a$epsilon,accepted=a$accepted,
        converged=a$converged,mgcv_conv=a$mgcv_conv,iterations=if(is.null(a$iterations))NA_integer_ else a$iterations,
        finite_fit=a$finite_fit,maximum_outer_gradient=if(is.null(a$maximum_outer_gradient))NA_real_ else a$maximum_outer_gradient,
        error=if(is.null(a$error))NA_character_ else a$error,warnings=paste(unlist(a$warnings),collapse=" | "),elapsed_seconds=a$elapsed_seconds))))
}
score_validation <- function(wide,models,levels_by_sex) {
  if(job$fold=="all")return(data.table())
  keep<-which(wide$fold_id==as.integer(job$fold));out<-list();cursor<-0L
  if(!length(keep))stop("No validation records for the requested fold.")
  for(sex in c("female","male"))for(first in seq.int(1L,length(keep),by=10000L)) {
    rows<-keep[seq.int(first,min(length(keep),first+9999L))]
    q<-wide[rows];y<-q[[paste0("y_",sex)]];mu<-rep(NA_real_,length(rows))
    known<-q$city%in%levels_by_sex[[sex]]
    if(!is.null(models[[sex]])&&any(known)) {
      q[,city_fac:=factor(city,levels=levels_by_sex[[sex]])]
      mu[known]<-predict(models[[sex]],as.data.frame(q[known]),type="response",discrete=FALSE)
    }
    # An all-zero training FE has a boundary mean of zero; report rather than silently drop it.
    train_present<-q$city%in%unique(wide$city[wide$fold_id!=as.integer(job$fold)])
    mu[!known&train_present]<-0
    valid<-is.finite(mu)&mu>=0;dev<-rep(NA_real_,length(rows));dev[valid]<-v2_poisson_deviance(y[valid],mu[valid])
    cursor<-cursor+1L
    out[[cursor]]<-data.table(city=q$city,y=y,mu=mu,dev=dev,known=known,valid=valid)[,
      .(n=.N,n_positive=sum(y>0),n_scored=sum(valid),n_unscored=sum(!valid),
        n_boundary=sum(!known&valid),observed_sum=sum(y),predicted_sum=sum(mu[valid]),
        poisson_deviance_sum=sum(dev[valid]),squared_error_sum=sum((y[valid]-mu[valid])^2),
        absolute_error_sum=sum(abs(y[valid]-mu[valid]))),by=city][,sex:=sex]
  }
  ans<-rbindlist(out)[,lapply(.SD,sum),by=.(city,sex),.SDcols=c("n","n_positive","n_scored","n_unscored",
    "n_boundary","observed_sum","predicted_sum","poisson_deviance_sum","squared_error_sum","absolute_error_sum")]
  ans[,`:=`(poisson_deviance_mean=poisson_deviance_sum/n,observed_mean=observed_sum/n,
    predicted_mean=predicted_sum/n,calibration_ratio=observed_sum/predicted_sum,
    score_status=ifelse(n_unscored>0,"unscored_rows_present",ifelse(is.infinite(poisson_deviance_sum),"infinite_boundary_loss","complete")),
    fold=job$fold,period=job$period,scope=job$scope,family=job$family,k=job$k,variant=job$variant,task_id=ID)]
  ans
}
run <- function() {
  status_write("running");log_stage("START",ID)
  sources<-v2_sources_for_variant(job)
  identity<-c(sources$identity,list(engine=v2_record(SCRIPT),helper=v2_record(file.path(dirname(SCRIPT),"v2_model_common.R")),
    variant_adapter=v2_record(file.path(dirname(SCRIPT),"v2_expanded_sources.R")),solver_repair=V2_SOLVER_REPAIR))
  if(job$op=="smoke") {
    wide<-smoke_data(job$k);blocks<-sort(unique(wide$spatial_block_2km));fold_identity<-NULL
    scaling<-v2_training_scaling(if(job$fold=="all")wide else wide[fold_id!=as.integer(job$fold)],
      allow_constant=job$architecture=="independent")
  }else {
    if(!length(sources$items)) {
      reason<-if(job$variant=="bridge")"excluded_by_frozen_V1_bridge_cohort"else"no_observed_input_rows"
      v2_csv_write(data.table(city=job$city,sex=c("female","male"),variant=job$variant,status=reason),file.path(OUT,"TRAIN_STATUS.csv"))
      status_write("not_estimable",list(reason=reason));return(invisible(NULL))
    }
    cached<-v2_raw_cached(sources,job,make_data,log_stage)
    wide<-cached$data;identity$raw_cache<-cached$identity;rm(cached)
    attached<-v2_attach_folds(wide,job$root,job$fold!="all"||job$rep>0L)
    wide<-attached$data;blocks<-attached$blocks;fold_identity<-attached$identity;rm(attached)
    tr<-if(job$fold=="all")rep.int(TRUE,nrow(wide))else wide$fold_id!=as.integer(job$fold)
    if(!any(tr)||!any(!tr)&&job$fold!="all")stop("Empty training/validation fold.")
    scaling<-if(job$fold=="all")sources$scaling else v2_training_scaling(wide[tr],
      allow_constant=job$architecture=="independent")
  }
  wide<-v2_scale(wide,scaling)
  # Training scaling and the immutable raw cache are already complete. All later
  # model/CV paths use z_* controls; releasing raw copies saves ~385 MiB at N=2.8m.
  wide[,c(V2_FIELDS,intersect("spatial_block_5km",names(wide))):=NULL]
  gc()
  wide[,case_weight:=v2_weights(wide,job$rep,blocks)]
  v2_csv_write(scaling,file.path(OUT,"SCALING.csv"))
  train_selector<-if(job$fold=="all")rep.int(TRUE,nrow(wide))else wide$fold_id!=as.integer(job$fold)
  uniqueL<-sort(unique(wide$z_road_mean[train_selector]))
  if(length(uniqueL)<job$k)stop("Requested k exceeds training unique L; no silent k reduction.")
  knots<-as.numeric(quantile(uniqueL,seq(0,1,length.out=job$k),names=FALSE,type=7))
  if(anyDuplicated(knots))stop("Duplicate training L knots.")
  source_rows<-nrow(wide)
  # Park the immutable scaled/fold/weight frame once. The optimizer must not keep
  # both the full paired wide table and a copied sex training table alive.
  work_frame<-file.path(job$root,"cache","active_wide",paste0(ID,"__pid",Sys.getpid(),".rds"))
  v2_rds_write(wide,work_frame)
  on.exit({
    if(file.exists(work_frame)&&!v2_retry(function()unlink(work_frame)==0L))
      log_stage("WORK_FRAME_CLEANUP_PENDING",work_frame)
  },add=TRUE)
  rm(train_selector,uniqueL)
  if(exists("tr",inherits=FALSE))rm(tr)
  gc()
  models<-summaries<-levels_by_sex<-controls_by_sex<-list();coefs<-scores<-statuses<-drops<-list()
  grid_residuals<-k_checks<-concurvities<-solver_attempts<-list()
  for(sex in c("female","male")) {
    if(!exists("wide",inherits=FALSE)){wide<-readRDS(work_frame);setDT(wide)}
    prepared<-list2env(prepare_sex(wide,sex),parent=emptyenv())
    rm(wide);gc()
    statuses[[sex]]<-copy(prepared$training)[,`:=`(sex=sex,period=job$period,task_id=ID,
      sex_model_status=prepared$status)]
    if(nrow(prepared$dropped))drops[[sex]]<-copy(prepared$dropped)[,sex:=sex]
    controls_by_sex[[sex]]<-prepared$controls;levels_by_sex[[sex]]<-prepared$city_levels
    if(prepared$status!="ready") {
      models[sex]<-list(NULL);summaries[[sex]]<-list(status=prepared$status)
      log_stage("NOT_ESTIMABLE",sex,prepared$status)
    }else {
      result<-fit_sex(prepared,sex,knots)
      models[[sex]]<-result$model;summaries[[sex]]<-result$summary;coefs[[sex]]<-result$coef;scores[[sex]]<-result$training_score
      solver_attempts[[sex]]<-result$solver_attempts
      if(nrow(result$grid_residuals))grid_residuals[[sex]]<-result$grid_residuals
      if(!is.null(result$diagnostics)) {
        k_checks[[sex]]<-result$diagnostics$k_check;concurvities[[sex]]<-result$diagnostics$concurvity
      }
      rm(result)
    }
    rm(prepared);gc()
    v2_csv_write(rbindlist(statuses,fill=TRUE),file.path(OUT,"TRAIN_STATUS.csv"))
  }
  metadata<-list(job=job,task_id=ID,scope=job$scope,architecture=job$architecture,city=job$city,
    period=job$period,family=job$family,k=job$k,fold=job$fold,rep=job$rep,
    phase=if(job$op=="smoke")"synthetic_smoke"else if(job$rep>0L)"paired_FRW_refit"else
      if(job$variant=="bridge")"same_target_bridge"else"expanded_observed_mean_analysis",
    variant=job$variant,cohort_variant=sources$cohort_variant,declared_cities=as.list(sources$declared),
    cohort_cities=as.list(sources$cohort),bridge_cities=if(job$variant=="bridge")as.list(sources$cohort)else NULL,
    city_levels_by_sex=levels_by_sex,control_columns_by_sex=controls_by_sex,L_knots=as.list(knots),
    reference=sources$reference,source_identity=identity,fold_map=fold_identity,
    scaling_policy=if(job$fold=="all")"Frozen V1 period scaling"else"Training fold only, unweighted mean/sample SD",
    zero_policy="Keep original finite encoded zeros; no NA-to-zero or new missing-row completion.",
    response="Raw nonnegative length-allocated archived PA mean, including encoded zero; not unique people.",
    independent_sex_fits="Sequential separate-sex fits implement the same sex-specific coefficients and smooths, with sex-city intercepts.",
    FRW=if(job$rep>0)list(seed=20261102L+job$rep,registry="FOLD_MAP global sorted city::2km blocks",
      distribution="Independent Gamma(1,1), normalized by global block mean; shared across sexes, periods, scopes and k.")else NULL,
    intervals="Model-based Vp in model; formal paired complete-refit uncertainty is separate.")
  metadata$boundary_cities_by_sex<-lapply(statuses,function(x)sort(unique(x$city[
    !is.na(x$training)&x$training&x$n_positive==0])))
  metadata$boundary_prediction_policy<-list(raw_mean="Zero only for explicitly registered all-zero training city-sex FE boundary cells",
    relative_curve="Undefined for a zero mean/reference; never fabricate a zero delta or unit ratio",
    other_unpredictable_states="Not zero; preserve missing/failed status separately")
  metadata$method_status<-list(fit=list(status=if(all(vapply(summaries,function(x)identical(x$status,"fitted"),logical(1))))"complete"else"partially_estimable",
    engine="Sequential sex-specific discrete bam; explicit G family; response-only prefit mf; compact prediction checked"),
    numerical_equivalence=list(status=if(job$op=="smoke")"per_sex_smoke_results"else"separate_synthetic_smoke_required",by_sex=lapply(summaries,`[[`,"numerical_qa")),
    grid_residuals=list(status=if(length(grid_residuals))"available_for_fitted_city_sex_cells_see_by_sex_status"else"not_available_or_not_requested",
      by_sex=lapply(summaries,`[[`,"grid_residuals"),
      omitted_all_zero_FE_cells="TRAIN_STATUS.csv records boundaries; no undefined Pearson residual fabricated"),
    model_diagnostics=list(status=if(length(k_checks)||length(concurvities))"fixed_seed_training_subsample"else"not_available_or_not_requested",
      by_sex=lapply(summaries,`[[`,"subsample_diagnostics")),
    spatial_test=list(status="awaiting_coordinate_join_and_grid_mean_spatial_diagnostics"),
    complete_refit_intervals=list(status=if(job$rep>0L)"this_paired_FRW_replicate_only"else"separate_stage_not_in_model_based_Vp"))
  if(job$op!="smoke") {
    omitted<-setdiff(if(job$architecture=="independent")job$city else sources$declared,sources$cohort)
    if(length(omitted))statuses$cohort_registry<-CJ(city=omitted,sex=c("female","male"))[
      ,`:=`(period=job$period,task_id=ID,
        status=if(job$variant=="bridge")"excluded_by_frozen_V1_bridge_cohort"else"no_observed_input_rows",training=NA,
        sex_model_status=if(job$variant=="bridge")"bridge_cohort_exclusion_not_invented_zero"else"empty_input_not_invented_zero")]
  }
  statuses<-lapply(statuses,function(x)x[,variant:=job$variant])
  v2_csv_write(rbindlist(statuses,fill=TRUE),file.path(OUT,"TRAIN_STATUS.csv"))
  bundle<-list(models=models,metadata=metadata,scaling=scaling,reference=sources$reference)
  v2_rds_write(bundle,file.path(OUT,"MODEL.rds"))
  v2_csv_write(rbindlist(coefs,fill=TRUE),file.path(OUT,"COEF.csv"))
  v2_csv_write(rbindlist(scores,fill=TRUE),file.path(OUT,"TRAIN_SCORE.csv"))
  v2_csv_write(rbindlist(drops,fill=TRUE),file.path(OUT,"REMOVED_CONTROLS.csv"))
  if(length(grid_residuals))v2_csv_write(rbindlist(grid_residuals,fill=TRUE),file.path(OUT,"GRID_RESIDUALS.csv"))
  if(length(k_checks))v2_csv_write(rbindlist(k_checks,fill=TRUE),file.path(OUT,"K_CHECK.csv"))
  if(length(concurvities))v2_csv_write(rbindlist(concurvities,fill=TRUE),file.path(OUT,"CONCURVITY.csv"))
  if(length(solver_attempts))v2_csv_write(rbindlist(solver_attempts,fill=TRUE),file.path(OUT,"SOLVER_ATTEMPTS.csv"))
  if(job$fold!="all") {
    wide<-readRDS(work_frame);setDT(wide)
    v2_csv_write(score_validation(wide,models,levels_by_sex),file.path(OUT,"CV_SCORE.csv"))
    rm(wide);gc()
  }
  complete<-all(vapply(summaries,function(x)identical(x$status,"fitted"),logical(1)))
  v2_json_write(list(status=if(complete)"complete"else"partially_estimable",task_id=ID,job=job,
    models=summaries,MODEL=v2_record(file.path(OUT,"MODEL.rds")),source_identity=identity,
    model_results_are_synthetic=job$op=="smoke",source_rows=source_rows,
    unknowns_filled_zero=FALSE,smearing_used=FALSE),file.path(OUT,"METRICS.json"))
  status_write(if(complete)"complete"else"partially_estimable");log_stage("DONE",ID)
}
tryCatch(run(),error=function(e) {
  status_write("failed",list(error=conditionMessage(e)));log_stage("FAILED",conditionMessage(e));stop(e)
})
