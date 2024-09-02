require(deSolve)
library(DescTools)

#Reflects 8/30 parameter name changes
SIR.joint <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    #Human population attributes
    N <- S+Ib+Ip+Ipb+Rp+Rpb
    Nb <- Ib+Ipb+Rpb
    Np <- Ip+Ipb
    
    #Mosquito population and related attributes
    M = m*(1-(1-alpha.m)*Nb)
    beta = (M*(a^2)*b*c*exp(-g*v))/g
    
    #Force of infection for behavior
    lambda.b <- phi.e*(phi.d*(Np)+phi.m*(M) + phi.s*(Nb))
    
    #Force of infection for pathogen without behavior
    lambda.p <- beta*Np
    
    #Force of infection for pathogen if doing behavior
    lambda.pb <- (alpha.b^2)*lambda.p
    
    dS <- omega*Ib + xi*Rp - lambda.p*S - lambda.b*S
    dIp <- lambda.p*S + omega*Ipb - lambda.b*Ip - gamma*Ip
    dIpb <- lambda.b*Ip + lambda.pb*Ib - gamma*Ipb - omega*Ipb
    dIb <- lambda.b*S - lambda.pb*Ib - omega*Ib + xi*Rpb
    dRp <- gamma*Ip - xi*Rp - lambda.b*Rp + omega*Rpb
    dRpb <- lambda.b*Rp + gamma*Ipb - xi*Rpb - omega*Rpb
    
    return(list(c(dS, dIp, dIpb, dIb, dRp, dRpb)))
  })
}


#initial values
initial_state <- c(S=.999, Ip=.001, Ipb=0, Ib=0, Rp=0, Rpb=0)

#5 year long simulation
times <- 0:(365*5)

phi.m = 0.01 #weight of mosquitoes on control participation
phi.d = 0.01 #weight of disease on control participation
phi.s = 0.01 #weight of social influence on control participation

alpha.b = 0.569 #1-0.431 #personal protection
alpha.m = .53 # remaining habitat
omega = 1/15 #duration of behavior
gamma = 1/7 #duration of dengue infection
xi = 1/180 #duration of dengue immunity

#mosquito parameters
m = 1 #ratio of mosquitoes to humans (baseline)
v = 1/14 #1/eip
g = 0.18 #mosquito mortality rate
a = 0.76 #mosquito biting rate
b = 0.3 #probability of mosquito to human infection
c = 0.3 #probability of human to mosquito infection


### To run a simulation without behavior turned on (phi.e = 0)
params = c(gamma,phi.m,phi.d,phi.s,alpha.b,omega,xi,m,v,g,a,b,c,alpha.m,phi.e=0)
model <- ode(initial_state, times, SIR.joint, params)

### To run a simulation with behavior turned on (phi.e = 1)
params.b = c(gamma,phi.m,phi.d,phi.s,alpha.b,omega,xi,m,v,g,a,b,c,alpha.m,phi.e=1)
model.b <- ode(initial_state, times, SIR.joint, params.b)

#Verify equilibrium value - behavior-free model
print(tail(model[,"Ip"]))

#Verify equilibrium value - model with behavior
print(tail(model.b[,"Ip"])+tail(model.b[,"Ipb"]))

#Equilibrium behavior prevalence - model with behavior - 0.140
print(tail(model.b[,"Ib"])+tail(model.b[,"Ipb"])+tail(model.b[,"Rpb"]))

#To calculate cumulative annual incidence for the final year of the simulation

#model without behavior
x.m <- tail(model[,"time"],365)
y.m <-tail(model[,"Ip"],365)
AUC(x.m,y.m)/1000 #to account for population of 1000

x.m.b <- tail(model.b[,"time"],365)
y.m.b <-tail(model.b[,"Ip"],365)+tail(model.b[,"Ipb"],365)
AUC(x.m.b,y.m.b)/1000 #to account for population of 1000

# To plot model simulations with and without dynamic behavior

par(mfrow=c(1,2))

par(mar = c(4.5,4,1,.5))
plot(model[,"Ip"], type = "l", lwd = 2, ylab = "Dengue prevalence", xlab = "Time (Days)",
     col = "darkgreen",lty = 2, ylim = c(0,0.17))
lines(model.b[,"Ip"]+model.b[,"Ipb"],lwd = 2, col = "black")
legend(615,0.15,legend = c("Without behavior", "With behavior"),col = c("darkgreen","black"), 
       lwd = 2, lty = c(2,1),cex = 0.65)

par(mar = c(4.5,4,1,.5))
plot(model.b[,"Ib"]+model.b[,"Ipb"]+model.b[,"Rpb"], col = "black", type = "l", lwd = 2,
     ylab = "Behavior prevalence", xlab = "Time (Days)",ylim = c(0,0.17))

plot(model.b[,"Ib"]+model.b[,"Ipb"]+model.b[,"Rpb"], col = "black", type = "l", lwd = 2,
     ylab = "Behavior prevalence", xlab = "Time (Days)")
