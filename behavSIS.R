#Script to explore solutions to SIS model and plot solution stability
library(ggplot2)
library(ggpubr)

#Expression for negative solution
sol1 <- function(omega,M,phi.m,alpha.m,phi.s){
  s1 = (-omega - 2*M*phi.m + M*alpha.m*phi.m + phi.s - 
    sqrt(-4*M*phi.m*(M*phi.m - M*alpha.m*phi.m - phi.s) + (-omega - 2*M*phi.m + M*alpha.m*phi.m + phi.s)^2))/
    (2*(-M*phi.m + M*alpha.m*phi.m + phi.s))
    return(s1)
}

#Expression for positive solution
sol2 <- function(omega,M,phi.m,alpha.m,phi.s){
  s2 = (-omega - 2*M*phi.m + M*alpha.m*phi.m + phi.s +
          sqrt(-4*M*phi.m*(M*phi.m - M*alpha.m*phi.m - phi.s) + (-omega - 2*M*phi.m + M*alpha.m*phi.m + phi.s)^2))/
    (2*(-M*phi.m + M*alpha.m*phi.m + phi.s))
  return(s2)
}

#Calculate first derivative
firstDeriv <- function(Ib,omega,M,phi.m,alpha.m,phi.s){
  d1 <- -omega - M*(1 - Ib*(1 - alpha.m))*phi.m - Ib*phi.s + 
    (1 - Ib)*(M*(-1 + alpha.m)*phi.m + phi.s)
  return(d1)
}

#Verify value of denominator in positive solution (to ensure real number)
checkDenom <- function(omega,M,phi.m,alpha.m,phi.s){
  d <- (M*phi.m*(1-alpha.m))/phi.s
  return(d)
}

#take the positive solution and do linear w ti

M = 1
n = 100 #For plotting purposes

alpha.m = 0.53
omega = 1/15 #Behavior recovery rate


gm = rep(seq(0,1, length=n), each=n)
gs = rep(seq(0,1, length=n),n)

#To generate positive and negative solutions, calculate derivative and denominator
df = data.frame(gm,gs)
df$s1 <- sol1(omega,M,phi.m=df$gm,alpha.m,phi.s=df$gs)
df$s2 <- sol2(omega,M,phi.m=df$gm,alpha.m,phi.s=df$gs)
df$d1 <- firstDeriv(Ib=df$s2,omega,M,phi.m=df$gm,alpha.m,phi.s=df$gs)
df$denom <- checkDenom(omega,M,phi.m=df$gm,alpha.m,phi.s=df$gs)

#Calculate positive solutions for various durations of behavior
df$o1 <- sol2(omega=1,M,phi.m=df$gm,alpha.m,phi.s=df$gs)
df$o2 <- sol2(omega=1/5,M,phi.m=df$gm,alpha.m,phi.s=df$gs)
df$o3 <- sol2(omega=1/30,M,phi.m=df$gm,alpha.m,phi.s=df$gs)

#Generate plots of behavior equilibria across changing duration
plot1<- ggplot(df, aes(x = gm, y = gs, fill = o1)) + scale_fill_viridis_b(limits=c(0.0,1.0))+ geom_tile() + theme_classic() +
  xlab(expression(alpha.m [M])) + ylab(expression(alpha.m [S])) + labs(fill = "") + ggtitle("\u03C9 = 1") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

plot2<- ggplot(df, aes(x = gm, y = gs, fill = o2)) + scale_fill_viridis_b(limits=c(0.0,1.0))+ geom_tile() + theme_classic() +
  xlab(expression(alpha.m [M])) + ylab(expression(alpha.m [S])) + labs(fill = "") + ggtitle("\u03C9 = 1/5") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

plot3 <-ggplot(df, aes(x = gm, y = gs, fill =o3)) + scale_fill_viridis_b(limits=c(0.0,1.0))+ geom_tile() + theme_classic() +
  xlab(expression(alpha.m [M])) + ylab(expression(alpha.m [S])) + labs(fill = "") + ggtitle("\u03C9 = 1/30") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

plot4 <- ggplot(df, aes(x = gm, y = gs, fill =o3)) + scale_fill_viridis_b(limits=c(0.0,1.0))+ geom_tile() +
  lims(x = c(0,0), y = c(0,0))+
  theme_void()+
  theme(legend.position = c(0.5,0.5), legend.title = element_blank())

ggarrange(plot1,plot2,plot3,plot4)