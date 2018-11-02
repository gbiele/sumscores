plot_densities = function(x, add = F) {
  dd = apply(x,2,function(x) density(x[!is.na(x)]))
  if (add == F) {
    ylim = range(sapply(dd, function(x) range(x$y)))
    xlim = range(sapply(dd, function(x) range(x$x)))
    plot(0,type = "n", ylim = ylim, xlim = xlim,ylab = "", xlab = ifelse(add,"b_hat","z"))
    legend("topleft",lty = 1, col = cols, legend = colnames(x),bty = "n")
  }
  x = sapply(1:ncol(x), function(x) lines(dd[[x]], col = cols[x], lty = ifelse(add,2,1)))
  return()
}


my_sim = function(mu_phi_b) {
  mu = mu_phi_b[1]
  phi = mu_phi_b[2]
  b = mu_phi_b[3]
  X = matrix(rnorm(N*K),ncol = K)
  p = inv.logit(mu + X %*% b)
  Y = rbetabinom.ab(N,size = 22, shape1 = p*phi, shape2 = (1-p)*phi)
  
  Ybin = Y > median(Y)
  bb_fit = BBreg(Y~X, m = 22)
  lm_fit = lm(Y~X)
  lr_fit = glm(Ybin~X, family = "binomial")
  
  return(c(summary(bb_fit)$coefficients[-1,3],
           summary(lm_fit)$coefficients[-1,3],
           summary(lr_fit)$coefficients[-1,3],
           diff(inv.logit(cumsum(bb_fit$coefficients))),
           lm_fit$coefficients[2]/22,
           diff(inv.logit(cumsum(lr_fit$coefficients)))
           ))
}

my_scatter = function (X,Y,N) {
  s = sample(500,N)
  
  x_scale = coef(lm(c(15,22)~range(X)))
  scaled_X = x_scale[1] + X*x_scale[2]
  usr = par("usr")
  yr <- (usr[4] - usr[3]) / 27
  ymax <- usr[4] - yr
  y_scale = coef(lm(ymax*c(.6,1)~range(Y)))
  scaled_Y = y_scale[1] + Y*y_scale[2]
  points(scaled_X[s],scaled_Y[s],pch = 16, col = adjustcolor("grey", alpha = .4))
  
  lmx = coef(lm(scaled_Y~scaled_X))
  lines(c(15,22),lmx[1] + lmx[2]*c(15,22), col = adjustcolor("red",.5), lwd = 3)
}

dist_plotter = function(mu,b,phi, xaxt = "s") {
  beta = (logit(inv.logit(mu)+b)-mu)
  X = matrix(rnorm(5000),ncol = K)
  p = inv.logit(mu + X %*% beta)
  
  sim_data = rbetabinom.ab(5000,
                           size = 22,
                           shape1 = p*phi,
                           shape2 = (1-p)*phi)
  hist(sim_data,
       breaks = (0:23)-.5,
       main = paste("b=",b,", mu=",inv.logit(mu),", phi=",phi),
       xlab = "sum score", ylab = "", yaxt = "n", xaxt = xaxt)
  my_scatter(X,sim_data,N)
  lines(11.5+c(0,-b*22),rep(par("usr")[4],2), lwd = 3, col = "red")
}