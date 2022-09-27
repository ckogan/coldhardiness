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
    # re_mat <- newdata[,re.form, drop = F] #%>% unique(); could do unique here for speed savings
    # unames <- lapply(1:nrow(re_mat), function(i) paste0(names(re.form), "[", re_mat[i,], "]"))
    #
    # u_samp <- lapply(unames, function(p) object$sflat[,dimnames(object$sflat)[[2]] %in% p])
    u_samp <- u_pars(object, newdata, re.form)
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
u_pars <- function(object, newdata, re.form) {
  re_mat <- newdata[,re.form, drop = F] #%>% unique(); could do unique here for speed savings
  unames <- lapply(1:nrow(re_mat), function(i) paste0(names(re.form), "[", re_mat[i,], "]"))
  lapply(unames, function(p) object$sflat[,dimnames(object$sflat)[[2]] %in% p])
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

#' @export
predict_marginal_sim <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL, vc.form.marginal = NULL, type=c("link", "response"), stat = NULL, parallel = T, B = 100, method = c("within", "side")) {
  if(parallel)
    `%fdo%` <- `%dopar%`
  else
    `%fdo%` <- `%do%`
  #lt <- vector("numeric", nrow(newdata))
  samples <- foreach(i = 1:nrow(newdata), .export = c("object","newdata", "re.form.conditional", "re.form.marginal"), .packages = c("coldhardiness")) %fdo% {
    cat(paste0(i, "/", nrow(newdata), "\n"))

    lp <- predict(object, newdata = newdata, re.form = re.form.conditional, type = "link", stat = NULL)
    if(!is.na(re.form.marginal))
      u_samp <- u_pars(object, newdata[i,], re.form.marginal)[[1]] #because only doing for one row
    pred_conditional <- array(NA, dim = c(dim(lp), B))
    for(j in 1:ncol(lp)) {
      vc <- object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal]^2
      re <- rnorm(B, 0, sqrt(rowSums(vc)))
      if(!is.na(re.form.marginal)) {
        re2 <- apply(u_samp, 2, function(x) sample(x, B, replace = T))
        re2 <- rowSums(re2)
        re <- re + re2
      }
      #lp_j <- sweep(re, 1, lp[,j], `+`)
      lp_j <- t(outer(re, lp[,j], FUN = `+`))
      if(type == "link") {
        pred <- lp_j
      } else if(type == "response") {
        pred <- 1 / (1 + exp(-(lp_j)))
      }
      pred_conditional[,j,] <- pred
    }
    apply(pred_conditional, c(1,2), stat)
  }

  do.call(partial(abind, along = 3), samples)
}

#' @export
predict_marginal_sim2 <- function(object, newdata = NULL, re.form.conditional = NULL, vc.form.marginal = NULL, type=c("link", "response")) {
    lp <- predict(object, newdata = newdata, re.form = re.form.conditional, type = "link", stat = NULL)
    vc <- object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal]^2
    re <- rnorm(ncol(lp), 0, sqrt(rowSums(vc)))
    lp <- sweep(lp, 2, re, `+`)
    if(type == "link") {
      pred <- lp
    } else if(type == "response") {
      pred <- 1 / (1 + exp(-(lp)))
    }
    newdata$pred <- rowMeans(pred)
    newdata
}

#' @export
augment_marginal_lt <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL, vc.form.marginal = NULL,  f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_names = "LT50", center = median, se_fit = NULL, std_ftemp_seq = seq(-2, 2, by = 0.005), parallel = T, B = 100) {
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
    if(!is.na(re.form.marginal))
      u_samp <- u_pars(object, newdata[i,], re.form.marginal)[[1]] #because only doing for one row
    p_conditional <- array(NA, dim = c(dim(lp), B))
    for(j in 1:ncol(lp)) {
      #re <- matrix(rnorm(nrow(lp)*B, 0, rowSums(object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal])), ncol = B)
      print(sqrt(rowSums(object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal]^2)))
      re <- rnorm(B, 0, sqrt(rowSums(object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal]^2)))
      if(!is.na(re.form.marginal)) {
        re2 <- apply(u_samp, 2, function(x) sample(x, B, replace = T))
        # print(apply(u_samp, 2, sd))
        # re2 <- sum(colMeans(u_samp))
        # re2 <- do.call(cbind, re2)
        # re2 <- apply(object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal], 2, function(x) sample(x, B, replace = T))
        re2 <- rowSums(re2)
        re <- re + re2
      }
      #lp_j <- sweep(re, 1, lp[,j], `+`)
      lp_j <- t(outer(re, lp[,j], FUN = `+`))
      p_conditional[,j,] <- 1 / (1 + exp(-(lp_j)))
    }
    p_marginal <- apply(p_conditional, c(1,2), mean)
    sapply(p, function(p_i) apply(p_marginal, 2, function(x) std_ftemp_seq[which(x > p_i)[1]])) # S x p
  }

  lt_samples <- do.call(partial(abind, along = 3), lt_samples)
  lt_samples <- f(lt_samples)
  lt <- apply(lt_samples, c(2,3), center)
  if(!is.null(se_fit)) {
    lt_se <- apply(lt_samples, c(2,3), se_fit)
  }

  for(i in 1:length(p)) {
    newdata[[paste0("LT", lt_names[i])]] <- lt[i,]
    if(!is.null(se_fit)) {
      newdata[[paste0("SE",lt_names[i])]] <- lt_se[i,]
    }
  }
  newdata
}


#' @export
augment_marginal_lt_simple <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL, vc.form.marginal = NULL,  f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_names = "LT50", std_ftemp_seq = seq(-2, 2, by = 0.005), parallel = T, B = 100) {
  # beta = NULL,
  # if(is.null(beta))
  if(parallel)
    `%fdo%` <- `%dopar%`
  else
    `%fdo%` <- `%do%`
  newdata[["std_ftemp"]] <- 0
  newdata <- unique(newdata)
  re <- rnorm(dim(object$sflat)[1], 0, sqrt(rowSums(object$sflat[,dimnames(object$sflat)[[2]] %in% vc.form.marginal]^2)))

  #lt <- vector("numeric", nrow(newdata))
  lt_samples <- foreach(i = 1:nrow(newdata), .export = c("object","newdata", "re.form.conditional", "re.form.marginal", "p", "std_ftemp_seq"), .packages = c("coldhardiness")) %fdo% {
    cat(paste0(i, "/", nrow(newdata), "\n"))
    data <- newdata[rep(i, length(std_ftemp_seq)),]
    data$std_ftemp <- std_ftemp_seq

    lp <- predict(object, newdata = data, re.form = re.form.conditional, type = "link", stat = NULL)
    if(!is.na(re.form.marginal))
      u_samp <- u_pars(object, newdata[i,], re.form.marginal)[[1]] #because only doing for one row

    lp <- sweep(lp, 2, re, `+`)
    p_conditional <- 1 / (1 + exp(-(lp)))

    p_marginal <- rowMeans(p_conditional)#apply(p_conditional, c(1,2), mean)
    sapply(p, function(p_i) std_ftemp_seq[which(p_marginal > p_i)[1]]) # S x p
  }

  lt_samples <- do.call(c, lt_samples)
  lt <- f(lt_samples)

  if(length(p) > 1) {
    for(i in 1:length(p)) {
      newdata[[paste0("LT", lt_names[i])]] <- lt[i,]
    }
  } else {
    newdata[[paste0("LT", lt_names)]] <- lt
  }

  newdata
}



#
# augment_marginal_lt <- function(object, newdata = NULL, re.form.conditional = NULL, re.form.marginal = NULL, vc.form.marginal = NULL,  f = compose(c_to_f, unstd_ftemp), p = 0.5, lt_names = "LT50", center = median, se_fit = NULL, std_ftemp_seq = seq(-2, 2, by = 0.005), parallel = T, B = 100) {
#   # newdata[["std_ftemp"]] <- 0
#   newdata <- unique(newdata)
#   lt <- lt_se <- list()
#   for(i in 1:nrow(newdata)) {
#     cat(paste0(i, "/", nrow(newdata)))
#     data <- newdata[rep(i, length(std_ftemp_seq)),]
#     data$std_ftemp <- std_ftemp_seq
#     # data <- newdata[rep(1:nrow(newdata), each = length(std_ftemp_seq)),]
#     # data$std_ftemp <- rep(std_ftemp_seq, nrow(newdata))
#     samples <- predict_marginal_sim(object, newdata = data, re.form.conditional = re.form.conditional, re.form.marginal = re.form.marginal, vc.form.marginal = vc.form.marginal, type="response")
#     lt_samples <- lapply(samples, function(p_marginal) sapply(p, function(p_i) apply(p_marginal, 2, function(x) std_ftemp_seq[which(x > p_i)[1]])))
#     lt_samples <- do.call(partial(abind, along = 3), lt_samples)
#     lt_samples <- f(lt_samples)
#     lt[[i]] <- apply(lt_samples, c(2,3), center)
#     if(!is.null(se_fit)) {
#       lt_se[[i]] <- apply(lt_samples, c(2,3), se_fit)
#     }
#   }
#   browser()
#
#
#   for(i in 1:length(p)) {
#     newdata[[paste0("LT", lt_names[i])]] <- lt[i,]
#     if(!is.null(se_fit)) {
#       newdata[[paste0("SE",lt_names[i])]] <- lt_se[i,]
#     }
#   }
#   newdata
# }

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

