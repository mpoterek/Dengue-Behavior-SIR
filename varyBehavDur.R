#Script to plot disease and behavior equilibrium values across omega parameter space
#where omega denotes the duration of behavior "infection"; we also evaluate this effect 
#while removing drivers of behavior one by one

source("jointModel.R") #To call model function

#Range of values for omega
omega.range = (seq(1,30,by=1))

#Presence of one phi at a time
all.eq = c()
gm.eq = c()
gs.eq=c()
gd.eq=c()

all.b = c()
gm.b = c()
gs.b = c()
gd.b = c()

for (i in 1:length(omega.range)){
  p = c(gamma,phi.m,phi.d,phi.s,alpha.b,omega=1/omega.range[i],xi,m,v,g,a,b,c,phi.e=1)
  p.gm = c(gamma,phi.m,phi.d=0,phi.s=0,alpha.b,omega=1/omega.range[i],xi,m,v,g,a,b,c,phi.e=1)
  p.gs = c(gamma,phi.m=0,phi.d=0,phi.s,alpha.b,omega=1/omega.range[i],xi,m,v,g,a,b,c,phi.e=1)
  p.gd = c(gamma,phi.m=0,phi.d,phi.s=0,alpha.b,omega=1/omega.range[i],xi,m,v,g,a,b,c,phi.e=1)
  
  all.mod = ode(initial_state, times, SIR.joint, p)
  all.eq[i] = tail(all.mod[,"Ip"],1)+tail(all.mod[,"Ipb"],1)
  all.b[i] = tail(all.mod[,"Ib"],1)+tail(all.mod[,"Ipb"],1)+tail(all.mod[,"Rpb"],1)
  
  gm.mod = ode(initial_state, times, SIR.joint, p.gm)
  gm.eq[i] = tail(gm.mod[,"Ip"],1)+tail(gm.mod[,"Ipb"],1)
  gm.b[i] = tail(gm.mod[,"Ib"],1)+tail(gm.mod[,"Ipb"],1)+tail(gm.mod[,"Rpb"],1)
  
  gs.mod = ode(initial_state, times, SIR.joint, p.gs)
  gs.eq[i] = tail(gs.mod[,"Ip"],1)+tail(gs.mod[,"Ipb"],1)
  gs.b[i] = tail(gs.mod[,"Ib"],1)+tail(gs.mod[,"Ipb"],1)+tail(gs.mod[,"Rpb"],1)
  
  gd.mod = ode(initial_state, times, SIR.joint, p.gd)
  gd.eq[i] = tail(gd.mod[,"Ip"],1)+tail(gd.mod[,"Ipb"],1)
  gd.b[i] = tail(gd.mod[,"Ib"],1)+tail(gd.mod[,"Ipb"],1)+tail(gd.mod[,"Rpb"],1)
}

plot(omega.range,all.eq, type = "l", lwd = 2, ylim = c(0.016,0.022), 
     ylab = "equilibrium disease prevalence")
lines(omega.range,gm.eq, lwd = 2, col = "blue")
lines(omega.range,gs.eq, lwd = 2, col = "red")
lines(omega.range, gd.eq, lwd = 2, col = "goldenrod")

plot(omega.range,all.b, type = "l", lwd = 2, ylim = c(0.0,0.27),
     ylab = "equilibrium behavior prevalence")
lines(omega.range,gm.b, lwd = 2, col = "blue")
lines(omega.range,gs.b, lwd = 2, col = "red")
lines(omega.range, gd.b, lwd = 2, col = "goldenrod")
