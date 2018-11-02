library(HRQoL)
library(parallel)
library(VGAM)
library(boot)
library(TeachingDemos)
library(ggplot2)

par (mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)

cl <- makeCluster(4)
clusterEvalQ(cl, library(VGAM))
clusterEvalQ(cl, library(HRQoL))
clusterEvalQ(cl, library(boot))

cols = c("green","red","blue")


source("utils.R")

par(mfrow = c(1,1))
dist_plotter(mu = logit(.25),b = .025,phi = 15)



N = 100
K = 1
bs = c(.02, .03, .04, .05)
mus = logit(c(.1,.25,.5))
phis = c(5, 15)

cex.main_orig = par("cex.main")
mar_orig = par("mar")

par(mfrow = c(4,6),cex.main = .75, mar = c(0,0,3,0))
for (b in bs) {
  for (mu in mus)  {
    for  (phi in phis) {
      dist_plotter(mu = mu, b = b, phi = phi, xaxt = "n")
    }
  }
}
par(cex.main = cex.main_orig)
par(mar = mar_orig)


clusterExport(cl, c("my_sim","N","K"))

sims = vector(length = length(mus)*length(bs)*length(phis), mode = "list")
j = 0
all_stats = c()
for (mu in mus) {
  for (phi in phis) {
    for (b in bs) {
      beta = (logit(inv.logit(mu)+b)-mu)
      pars = matrix(rep(c(mu,phi,beta),2000),nrow = 3)
      stats = parApply(cl,pars,2,my_sim)
      stats = t(stats)
      colnames(stats) = c("beta-binomial","linear","logistic",
                          "me_bb","me_lin","me_lr")
      
      par(mfrow = c(2,2))
      p = inv.logit(mu + matrix(rnorm(N*K*10),ncol = K) %*% beta)
      hist(rbetabinom.ab(1000,
                         size = 22,
                         shape1 = p*phi,
                         shape2 = (1-p)*phi),
           main = "distribution of sum scores",
           xlab = "sum score", ylab = "", yaxt = "n")
      mtext(paste("mu = ",inv.logit(mu),", phi = ",phi,", b = ",b))
      plot_densities(stats[,1:3])
      title("distribution of z values")
      abline(v = 1.96, lty = 2, col = "grey")
      plot_densities(stats[,4:6])
      tmp = stats[,4:6]
      for (k in 1:3) tmp[stats[,k]<1.96,k] = NA
      plot_densities(tmp, add = T)
      title("estimated effects")
      legend("topright", col = "grey", legend = "only sign.")
      abline(v = b, lty = 2, col = "grey")
      barplot(colMeans(stats[,1:3] > 1.96),
              ylim = c(0,1), 
              col = adjustcolor(cols,alpha = .75),
              main = "power",
              xlab = "analysis model")
      
      all_stats = rbind(all_stats,
                        cbind(rep(mu,3),
                              rep(phi,3),
                              rep(b,3),
                              colMeans(stats[,1:3]),
                              colMeans(stats[,1:3] > 1.96),
                              colMeans(stats[,4:6])))
      colnames(all_stats) = c("mu","phi","b","z","power","b_hat")
      j = j+1
      sims[[j]] = list(mu = mu, b = b, phi = phi, stats = stats)
    }
  }
}

save(sims,all_stats,bs,mus,phis,file = "sims.Rdata")

stopCluster(cl)

df = data.frame(all_stats)
df$model = rep(c("bbr","l","lr"),length(bs)*3)
df$mu = inv.logit(df$mu)
ggplot(df,aes(x = b, y = power, col = model)) + geom_line() + facet_wrap(~phi+mu)
ggplot(df,aes(x = b, y = b_hat, col = model)) + geom_abline(intercept = 0, slope = 1, aes(colour='white')) + geom_line() + facet_wrap(~phi+mu)


ggplot(df,aes(x = b, y = power, col = model)) + geom_line() + facet_wrap(~phi+mu)
ggplot(df,aes(x = b, y = z, col = model)) + geom_line() + facet_wrap(~phi+mu)