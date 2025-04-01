####
# 

simTimeTrendsTreshold <- function(N,h2a,h2b,rg,r2a,r2b,prevs) {
covmatG=matrix(c(h2a,sqrt(h2a)*sqrt(h2b)*rg,sqrt(h2a)*sqrt(h2b)*rg,h2a),2)
liabG <- MASS::mvrnorm(N,c(0,0),Sigma = covmatG)
covmatE=matrix(c(1-h2a,0,0,1-h2b),2)
liabE <- MASS::mvrnorm(N,c(0,0),Sigma = covmatE)

liabTotal = liabG + liabE

# r2 = h2^2/(h2+varNoise)
# varNoise = h2^2/r2 - h2

covmatNoise=matrix(c(h2a^2/r2a - h2a,0,0,h2b^2/r2b - h2b),2)
Noise <- MASS::mvrnorm(N,c(0,0),Sigma = covmatNoise)
PGS= liabG+ Noise

PGS_scaled=apply(PGS, 2, scale)

df = data.table(PGS=PGS,PGS_scaled=PGS_scaled,l=liabTotal)

out = Reduce('rbind',lapply(prevs,function(prev){
df[l.V1>qnorm(1-prev),.(h2a,h2b,r2a,r2b,rg,prev,m_PGS_a=mean(PGS_scaled.V1),m_PGS_b=mean(PGS_scaled.V2))]}))
out}

df_out <- Reduce('rbind',lapply(c(0,.25,.5,.75), function(rg) 
  simTimeTrendsTreshold(N=5e5,h2a = .8,h2b = .8,rg = rg,r2a = .1,r2b = .1,prevs = seq(0.02,0.04,.005))))

g1 <- ggplot(data=melt(df_out,id.vars = c( "h2a", "h2b", "r2a", "r2b",   "rg",  "prev")),
             aes(x=prev,y=value,color=ifelse(variable=="m_PGS_a","Disease A","Disease B")))+
  geom_point()+geom_line()+facet_wrap(~paste("rg=",rg),ncol=4)+
  labs(title="h2a = h2b = 0.8, r2a = r2b = 0.1",y="PGS mean in cases of Disease A",color="PGS")+
    theme_bw()+scale_color_manual(values = c("blue","red"))+
  theme(plot.title = element_text(hjust = 0.5))


simTimeTrendsShifting <- function(N,h2a,h2b,rg,r2a,r2b,prevA=0.02,prevB=0.03,props_B_in_A) {
  covmatG=matrix(c(h2a,sqrt(h2a)*sqrt(h2b)*rg,sqrt(h2a)*sqrt(h2b)*rg,h2a),2)
  liabG <- MASS::mvrnorm(N,c(0,0),Sigma = covmatG)
  covmatE=matrix(c(1-h2a,0,0,1-h2b),2)
  liabE <- MASS::mvrnorm(N,c(0,0),Sigma = covmatE)
  
  liabTotal = liabG + liabE
  
  # r2 = h2^2/(h2+varNoise)
  # varNoise = h2^2/r2 - h2
  
  covmatNoise=matrix(c(h2a^2/r2a - h2a,0,0,h2b^2/r2b - h2b),2)
  Noise <- MASS::mvrnorm(N,c(0,0),Sigma = covmatNoise)
  PGS= liabG+ Noise
  
  PGS_scaled=apply(PGS, 2, scale)
  
   
  out = Reduce('rbind',lapply(props_B_in_A,function(prop_B_in_A){
    df = data.table(PGS=PGS,PGS_scaled=PGS_scaled,l=liabTotal)
    df[,caseA:=l.V1>qnorm(1-prevA)]
    df[caseA==T,ShiftedCaseA:=1]
    df[,caseB:=l.V2>qnorm(1-prevB)]
    # prop_B_in_A = sum(B!A)/(sum(A)+sum(B!A))
    # (sum(A)*prop_B_in_A+sum(B!A)*prop_B_in_A = sum(B!A)
    # (sum(A)*prop_B_in_A/ (1-*prop_B_in_A) = sum(B!A)
    
    N_B_extra= df[,sum(caseA)*prop_B_in_A/(1-prop_B_in_A)]
    N_B= df[,sum((caseA==0&caseB==1))]
    
    df[(caseA==0&caseB==1)&runif(N_B)<N_B_extra/N_B,ShiftedCaseA:=1]
    
    df[,.N,.(caseA,ShiftedCaseA,caseB)]
      df[ShiftedCaseA==1,.(h2a,h2b,r2a,r2b,rg,prevA,prevB,prop_B_in_A,m_PGS_a=mean(PGS_scaled.V1),m_PGS_b=mean(PGS_scaled.V2))]}))
  out}

df_out <- Reduce('rbind',lapply(c(0,.25,.5,.75), function(rg) 
  simTimeTrendsShifting(N=5e5,h2a = .8,h2b = .8,rg = rg,r2a = .1,r2b = .1,prevA = 0.02,prevB = 0.03,props_B_in_A = c(0,.125,.25,.375,.5) )))

g2 <- ggplot(data=melt(df_out,id.vars = c( "h2a", "h2b", "r2a", "r2b",   "rg",  "prop_B_in_A","prevA","prevB")),
             aes(x=prop_B_in_A,y=value,color=ifelse(variable=="m_PGS_a","Disease A","Disease B")))+
  geom_point()+geom_line()+facet_wrap(~paste("rg=",rg),ncol=4)+
  labs(title="h2a = h2b = 0.8, r2a = r2b = 0.1",y="PGS mean in cases of Disease A",color="PGS")+
  theme_bw()+scale_color_manual(values = c("blue","red"))+
  theme(plot.title = element_text(hjust = 0.5))



simTimeTrendsDetection <- function(N,h2a,h2b,rg,r2a,r2b,prevA=0.04,props_missed) {
  covmatG=matrix(c(h2a,sqrt(h2a)*sqrt(h2b)*rg,sqrt(h2a)*sqrt(h2b)*rg,h2a),2)
  liabG <- MASS::mvrnorm(N,c(0,0),Sigma = covmatG)
  covmatE=matrix(c(1-h2a,0,0,1-h2b),2)
  liabE <- MASS::mvrnorm(N,c(0,0),Sigma = covmatE)
  
  liabTotal = liabG + liabE
  
  # r2 = h2^2/(h2+varNoise)
  # varNoise = h2^2/r2 - h2
  
  covmatNoise=matrix(c(h2a^2/r2a - h2a,0,0,h2b^2/r2b - h2b),2)
  Noise <- MASS::mvrnorm(N,c(0,0),Sigma = covmatNoise)
  PGS= liabG+ Noise
  
  PGS_scaled=apply(PGS, 2, scale)
  
  
  out = Reduce('rbind',lapply(props_missed,function(prop_missed){
    df = data.table(PGS=PGS,PGS_scaled=PGS_scaled,l=liabTotal)
    df[,caseA:=l.V1>qnorm(1-prevA)]
    df[,ShiftedCaseA:=0]
    
    df[caseA==T,ShiftedCaseA:=1]

    
    N_A= df[,sum(caseA==1)]
    
    df[caseA==1&runif(N)<prop_missed,ShiftedCaseA:=0]
    
    df[,.N,.(caseA,ShiftedCaseA)]
    df[,.(h2a,h2b,r2a,r2b,rg,prevA,prop_missed,m_PGS_a=mean(PGS_scaled.V1),m_PGS_b=mean(PGS_scaled.V2)),ShiftedCaseA]}))
  out}

df_out <- Reduce('rbind',lapply(c(0,.25,.5,.75), function(rg) 
  simTimeTrendsDetection(N=5e6,h2a = .8,h2b = .8,rg = rg,r2a = .1,r2b = .1,prevA = 0.04,props_missed = c(0,.125,.25,.375,.5) )))

g3 <- ggplot(data=melt(df_out,id.vars = c( "h2a", "h2b", "r2a", "r2b",   "rg",  "prop_missed","prevA","ShiftedCaseA")),
             aes(x=prop_missed,y=value,color=ifelse(variable=="m_PGS_a","Disease A","Disease B"),linetype=ShiftedCaseA==T))+
  geom_point()+geom_line()+facet_wrap(~paste("rg=",rg),ncol=4)+
  labs(title="h2a = h2b = 0.8, r2a = r2b = 0.1",y="PGS mean in cases of Disease A",color="PGS",linetype="Case")+
  theme_bw()+scale_color_manual(values = c("blue","red"))+
  theme(plot.title = element_text(hjust = 0.5))+scale_x_reverse()
g3


