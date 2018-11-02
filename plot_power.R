library(boot)
library(VGAM)
library(TeachingDemos)
load("sims.Rdata")
K = 1
N = 100
breaks = -.5:22.5

df = data.frame(all_stats,row.names = NULL)
df$model = rep(c("bbr","l","lr"),length(bs)*3)
df$mu = inv.logit(df$mu)

cols = list(bbr = "green", l = "red", lr = "blue")
par(mfrow = c(2,3))
for (phi in phis)  {
  for   (mu in mus) {
    
    
    plot(0,type = "n", ylim = c(0,1), xlim = range(bs),
         ylab = "power",
         xlab = "b", bty = "n",xaxt = "n",
         main = paste0("mu=", inv.logit(mu), ", phi=",phi))
    axis(1,pos = 0)
    
    for (m in unique(df$model)) {
      tmp_stats = df[df$mu == inv.logit(mu) & df$phi == phi & df$model == m,]
      lines(tmp_stats$b, tmp_stats$power, col = cols[[m]])
    }
    
    y = rbetabinom.ab(5000,
                      size = 22,
                      shape1 = inv.logit(mu)*phi,
                      shape2 = (1-inv.logit(mu))*phi)
    h = hist(y,plot = F, breaks = breaks)
    
    x_scale = coef(lm(c(.04,.05)~range(h$mids)))
    scaled_X = x_scale[1] + h$mids*x_scale[2]
    y_scale = coef(lm(c(0,.33)~c(0,max(h$counts))))
    scaled_Y = y_scale[1] + h$counts*y_scale[2]
    
    lines(c(.0395,scaled_X),c(0,scaled_Y),"S", col = "darkred")
    legend("topleft", lty = 1, bty = "n",
           col = c("green","red","blue"),
           legend = c("beta binomial","linear","logistic"))
    
  }
}


par(mfrow = c(2,3))
for (phi in phis)  {
  for   (mu in mus) {
    
    
    plot(0,type = "n", ylim = range(bs), xlim = range(bs),
         ylab = "b-est",
         xlab = "b", bty = "n",xaxt = "n",
         main = paste0("mu=", inv.logit(mu), ", phi=",phi))
    axis(1,pos = 0.018)
    abline(0,1, col = "grey", lty = 2)
    for (m in unique(df$model)[1:2]) {
      tmp_stats = df[df$mu == inv.logit(mu) & df$phi == phi & df$model == m,]
      lines(tmp_stats$b, tmp_stats$b_hat, col = cols[[m]])
    }
    
    y = rbetabinom.ab(5000,
                      size = 22,
                      shape1 = inv.logit(mu)*phi,
                      shape2 = (1-inv.logit(mu))*phi)
    h = hist(y,plot = F, breaks = breaks)
    
    x_scale = coef(lm(c(.04,.05)~range(h$mids)))
    scaled_X = x_scale[1] + h$mids*x_scale[2]
    y_scale = coef(lm(c(.019,.03)~c(0,max(h$counts))))
    scaled_Y = y_scale[1] + h$counts*y_scale[2]
    
    lines(c(.0395,scaled_X),c(.018,scaled_Y),"S", col = "darkred")
    legend("topleft", lty = 1, bty = "n",
           col = c("green","red","blue"),
           legend = c("beta binomial","linear","logistic"))
    
  }
}

library(data.table)
dt = data.table(df)
for (b in unique(df$b)) {
  for (m in unique(df$model)[1:2]) {
    dt[model == m & b == b & z > 1.96]
  }
}

all_sims = c()
for (s in 1:length(sims)) {
  stats = data.table(sims[[s]]$stats)
  stats[,mu := round(inv.logit(sims[[s]]$mu),digits = 2)]
  stats[,b := sims[[s]]$b]
  stats[,phi := sims[[s]]$phi]
  all_sims = rbind(all_sims,stats)
}
all_sims$b = round(all_sims$b,digits = 2)

for (mb in c(0.02,0.03,0.04,0.05)) {
  d1f = density(all_sims[b == mb, me_bb])
  d2f = density(all_sims[b == mb, me_lin])
  d1 = density(all_sims[`beta-binomial` > 1.96 & b == mb, me_bb])
  d2 = density(all_sims[linear > 1.96 & b == mb, me_lin])
  plot(0,type = "n",
       xlim = c(-.05,.15),
       ylim = range(c(d1$y,d2$y,d1f$y,d2f$y)),
       xlab = "estimated effect b",
       ylab = "density",
       bty = "n",
       main = paste0("true b = ",mb))
  lines(d1, col = "green")
  lines(d2, col = "red")
  lines(d1f, col = "green", lty = 2)
  lines(d2f, col = "red", lty = 2)
  legend("topleft",
         col = c("green","red","grey","grey"),
         lty = c(1,1,1,2),
         bty = "n",
         legend = c("beta binomial","linear","all estimates", "significant estimates"))
  abline(v = mb)
}