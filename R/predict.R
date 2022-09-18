#' @export
chbmod <- function(object, rhsform) {
  sflat = flatten_stan_array(as.array(object))

  structure(list(model_pars = dimnames(sflat)[[2]],
                 rhsform = rhsform,
                 sflat = sflat),
            class = "chbmod")
}

#' @export
predict.chbmod <- function(object, newdata = NA, re.form=NULL, type=c("link", "response"), stat = NULL) {
  #re_pars = c(r_seasonXfield = "seasonXfield_ID", r_fieldXdate = "fieldXdate_ID"),
  X <- model.matrix(object$rhsform, newdata)
  lp <- X %*% beta_pars(object)

  if(!is.na(re.form)) {
    re_mat <- newdata[,re.form, drop = F] #%>% unique(); could do unique here for speed savings
    unames <- lapply(1:nrow(re_mat), function(i) paste0(names(re.form), "[", re_mat[i,], "]"))

    u_samp <- lapply(unames, function(p) object$sflat[,dimnames(object$sflat)[[2]] %in% p])
    zu_samp <- lapply(u_samp, rowSums)
    zu <- do.call(rbind, zu_samp)
    lp <- lp + zu
  }
  if(type == "link") {
    pred <- lp
  } else if(type == "response") {
    pred <- binomial()$linkinv(lp)
  }
  if(!is.null(stat))
    pred <- apply(pred, 1, stat)
  pred
}

#' @export
beta_pars <- function(obj) {
  t(samp(obj$sflat, "beta", partial = T))
}

#' @export
augment_lt <- function(object, newdata = weather_daily_default_fieldseason(), re.form=NULL, f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_names = "LT50", center = median, se_fit = NULL, xcol = 2) {
  # beta = NULL,
  # if(is.null(beta))
  newdata[["std_ftemp"]] <- 0
  newdata <- unique(newdata)
  lp <- predict(object, newdata = newdata, re.form = re.form, type = "link", stat = NULL)


  for(i in 1:length(p)) {
    lt_samples <- sweep(log(p[i]/(1-p[i])) - lp, 2, beta_pars(object)[xcol,],`/`)
    lt_samples <- f(lt_samples)
    lt <- apply(lt_samples, 1, center)
    newdata[[lt_names[i]]] <- lt
    if(!is.null(se_fit)) {
      newdata[[paste0(lt_names[i], "_se")]] <- apply(lt_samples, 1, se_fit)
    }
  }
  newdata
}

#' #' @export
#' augment_marginal_lt <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL,  f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_name = "LT50", center = median, se_fit = NULL, std_ftemp_seq = seq(-2, 2, by = 0.005)) {
#'   # beta = NULL,
#'   # if(is.null(beta))
#'   newdata[["std_ftemp"]] <- 0
#'   newdata <- unique(newdata)
#'   lt <- vector("numeric", nrow(newdata))
#'   for(i in 1:nrow(newdata)) {
#'     cat(paste0(i, "/", nrow(newdata), "\n"))
#'     data <- newdata[rep(i, length(std_ftemp_seq)),]
#'     data$std_ftemp <- std_ftemp_seq
#'
#'     lp <- predict(object, newdata = data, re.form = re.form.conditional, type = "link", stat = NULL)
#'
#'     re <- rnorm(ncol(lp), 0, rowSums(object$sflat[,dimnames(object$sflat)[[2]] %in% re.form.marginal]))
#'     lp <- sweep(lp, 2, re, `+`)
#'     p_conditional <- 1 / (1 + exp(-lp))
#'     p_marginal <- rowMeans(p_conditional)
#'     lt[i] <- std_ftemp_seq[which(p_marginal > p)[1]]
#'   }
#'   newdata[[lt_name]] <- f(lt)
#'   newdata
#' }


#' augment_marginal_lt <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL,  f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_name = "LT50", center = median, se_fit = NULL, std_ftemp_seq = seq(-2, 2, by = 0.005), parallel = T) {
#'   # beta = NULL,
#'   # if(is.null(beta))
#'   if(parallel)
#'     `%fdo%` <- `%dopar%`
#'   else
#'     `%fdo%` <- `%do%`
#'   newdata[["std_ftemp"]] <- 0
#'   newdata <- unique(newdata)
#'   #lt <- vector("numeric", nrow(newdata))
#'   lt <- foreach(i = 1:nrow(newdata), .export = c("object","newdata", "re.form.conditional", "re.form.marginal", "p", "std_ftemp_seq")) %fdo% {
#'     cat(paste0(i, "/", nrow(newdata), "\n"))
#'     data <- newdata[rep(i, length(std_ftemp_seq)),]
#'     data$std_ftemp <- std_ftemp_seq
#'
#'     lp <- predict(object, newdata = data, re.form = re.form.conditional, type = "link", stat = NULL)
#'
#'     re <- rnorm(ncol(lp), 0, rowSums(object$sflat[,dimnames(object$sflat)[[2]] %in% re.form.marginal]))
#'     lp <- sweep(lp, 2, re, `+`)
#'     p_conditional <- 1 / (1 + exp(-lp))
#'     p_marginal <- rowMeans(p_conditional)
#'     std_ftemp_seq[which(p_marginal > p)[1]]
#'   }
#'   lt <- do.call(c, lt)
#'   newdata[[lt_name]] <- f(lt)
#'   newdata
#' }

#' @export
augment_marginal_lt <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL,  f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_names = "LT50", center = median, se_fit = NULL, std_ftemp_seq = seq(-2, 2, by = 0.005), parallel = T, B = 100) {
  # beta = NULL,
  # if(is.null(beta))
  if(parallel)
    `%fdo%` <- `%dopar%`
  else
    `%fdo%` <- `%do%`
  newdata[["std_ftemp"]] <- 0
  newdata <- unique(newdata)
  #lt <- vector("numeric", nrow(newdata))
  lt_samples <- foreach(i = 1:nrow(newdata), .export = c("object","newdata", "re.form.conditional", "re.form.marginal", "p", "std_ftemp_seq"), .packages = c("coldhardiness")) %fdo% {
    cat(paste0(i, "/", nrow(newdata), "\n"))
    data <- newdata[rep(i, length(std_ftemp_seq)),]
    data$std_ftemp <- std_ftemp_seq

    lp <- predict(object, newdata = data, re.form = re.form.conditional, type = "link", stat = NULL)
    p_conditional <- array(NA, dim = c(dim(lp), B))
    for(j in 1:ncol(lp)) {
      re <- matrix(rnorm(nrow(lp)*B, 0, rowSums(object$sflat[,dimnames(object$sflat)[[2]] %in% re.form.marginal])), ncol = B)
      lp_j <- sweep(re, 1, lp[,j], `+`)
      p_conditional[,j,] <- 1 / (1 + exp(-(lp_j)))
    }
    p_marginal <- apply(p_conditional, c(1,2), mean)
    sapply(p, function(p_i) apply(p_marginal, 2, function(x) std_ftemp_seq[which(x > p_i)[1]])) # S x p
  }
  lt_samples <- do.call(partial(abind, along = 3), lt_samples)
  lt_samples <- f(lt_samples)
  lt <- apply(lt_samples, c(2,3), center)
  if(!is.null(se_fit))
    lt_se <- apply(lt_samples, c(2,3), se_fit)
  for(i in 1:length(p)) {
    newdata[[paste0("LT", lt_names[i])]] <- lt[i,]
    if(!is.null(se_fit)) {
      newdata[[paste0("SE",lt_names[i])]] <- lt_se[i,]
    }
  }
  newdata
}

#' @export
flatten_stan_array <- function(x) {
  lst <- lapply(1:dim(x)[2], function(i) x[,i,, drop = T])
  do.call(rbind, lst)
}

#' @export
flatten_stanfit <- function(obj) {
  if(class(obj) == "stanfit") {
    flatten_stan_array(as.array(obj))
  } else if(length(dim(obj)) == 3) {
    flatten_stan_array(obj)
  } else {
    obj
  }
}

# beta_pars <- function(obj) {
#   t(samp(obj$sflat, "beta", partial = T))
# }
#
# u_pars <- function(obj, ids) {
#   nm <- names(ids)
#   pars_lst <- lapply(1:length(ids), function(i) paste0(nm[i], "[", ids[i], "]"))
#   pars <- do.call(c, pars_lst)
#   obj$sflat[,obj$model_pars %in% pars]
# }


# predict.chbmod <- function(object, newdata = NA, re.form=NULL, type=c("link", "response", "lt"), stat = NULL, lt_par = 2, fresponse = NULL, p = 0.5) {
#   #re_pars = c(r_seasonXfield = "seasonXfield_ID", r_fieldXdate = "fieldXdate_ID"),
#   browser()
#   if(type == "lt"){
#     newdata[["std_ftemp"]] <- 0
#     newdata <- unique(newdata)
#   }
#
#   X <- model.matrix(object$rhsform, newdata)
#   lp <- X %*% beta_pars(object)
#
#   if(!is.na(re.form)) {
#     browser()
#     re_mat <- newdata[,re.form, drop = F] #%>% unique(); could do unique here for speed savings
#     unames <- lapply(1:nrow(re_mat), function(i) paste0(names(re.form), "[", re_mat[i,], "]"))
#
#     u_samp <- lapply(unames, function(p) object$sflat[,dimnames(object$sflat)[[2]] %in% p])
#     zu_samp <- lapply(u_samp, rowSums)
#     zu <- do.call(rbind, zu_samp)
#     lp <- lp + zu
#   }
#   browser()
#   if(type == "lt") {
#     pred <- sweep(log(p/(1-p)) - lp, 2, beta_pars(object)[lt_par,],`/`)
#   }
#   if(type == "response") {
#     pred <- binomial()$linkinv(pred)
#     if(!is.null(unstd))
#       pred <- unstd(pred)
#   }
#   if(!is.null(stat))
#     pred <- apply(pred, 1, stat)
#   browser()
#   pred
# }



# augment.chbmod <- function(object, newdata, ...) {
#   newdata[[".fitted"]] <- predict(object, newdata, ...)
#   newdata
# }
#
# betahat <- function(stanfit, stat = median) {
#   sflat <- flatten_stanfit(stanfit)
#   sbeta <- samp(sflat, "beta", partial = T)
#   beta <- apply(sbeta, 2, stat)
#   beta
# }
#
# sample_u <- function(stanfit, ids) {
#   sflat <- flatten_stanfit(stanfit)
#   nm <- names(ids)
#   pars_lst <- lapply(1:length(ids), function(i) paste0(nm[i], "[", ids[i], "]"))
#   pars <- do.call(c, pars_lst)
#   u <- sflat[,dimnames(sflat)[[2]] %in% pars]
#   u
# }
#
# sample_zu <- function(stanfit, ids) {
#   u <- sample_u(stanfit, ids)
#   rowMeans(u)
# }
#
# sample_lp <- function(stanfit, data, form, ids) {
#   sflat <- flatten_stanfit(stanfit)
#   sbeta <- samp(sflat, "beta", partial = T)
#   szu <- sample_zu(stanfit, ids)
#   xb <- lp_xb(data, form, sbeta)
#   sweep(xb, 2, szu, FUN = `+`)
# }
#
# doit <- function(stanfit, data, form, ids) {
#   sflat <- flatten_stanfit(stanfit)
#   lp <- sample_lp(stanfit, data, form, ids)
#   sbeta1 <- samp(sflat, "beta", partial = T)[,2]
#   LT(lp, sbeta1, 0.5)
# }
#
# uhat <- function(stanfit, ids, stat = median) {
#   sflat <- flatten_stanfit(stanfit)
#   nm <- names(ids)
#   pars_lst <- lapply(1:length(ids), function(i) paste0(nm[i], "[", ids[i], "]"))
#   pars <- do.call(c, pars_lst)
#   u <- apply(sflat[,dimnames(sflat)[[2]] %in% pars], 2, stat)
#   u <- as.list(u)
#   names(u) <- nm
#   u
# }
#
# lp_xb <- function(data, form, beta, xcol = 2) {
#   data <- data %>% mutate(std_ftemp = 0)
#   X <- model.matrix(form, data)
#   X[,xcol] <- 0
#   X %*% t(beta)
# }
#
#
# buhat <- function(stanfit, ids, stat = median) {
#   list(betahat = betahat(stanfit, stat = stat), uhat = uhat(stanfit, ids, stat = median))
# }
#
# lp_cond_int <- function(data, form, buhat, xcol = 2) {
#   xbint <- lp_xb_int(data, form, buhat$beta, xcol = xcol)
#   uint <- do.call(sum, buhat$uhat)
#   as.numeric(xbint + uint)
# }
#
# lp_xb_int <- function(data, form, beta, xcol = 2) {
#   data <- data %>% mutate(std_ftemp = 0)
#   X <- model.matrix(form, data)
#   X[,xcol] <- 0
#   as.numeric(X %*% beta)
# }
