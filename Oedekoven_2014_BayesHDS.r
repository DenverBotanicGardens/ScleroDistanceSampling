################################################################################################################################################

                                 # Author: Cornelia Oedekoven 
                                 
                                 # Centre for Research into Ecological and Environmental Modelling
                                 # University of St Andrews, St Andrews, Scotland
                                 # cornelia@mcs.st-and.ac.uk 
                                 
                                 # R code for RJMCMC algorithm for analysing covey data (repeated point transects with exact distance data)
                                 
                                 # Additional functions are given to adapt analysis to line transect and/or interval distance data 
                                 

################################################################################################################################################

# Data format required:

# for the detection function model L_y(\bmath{\theta}) (eqn (3))
# matrix covey.d: a n*6 matrix  (n = total number of detections, 6 columns: id, distance, site, year, type, state)
covey.d<-covey2[which(is.na(covey2$Distance)==F),]
length.d<-length(covey.d$Distance)

# for the density model L_n(\bmath{\beta}|\bmath{\theta}) (eqn (6))
# matrix n.jpr: a j*max(Tj) matrix with n_jpr (counts at site j, point p, visit r (Tj is a vector of numbers of measurements in the jth group)
# matrices with covariates: j*max(Tj) matrix with observed covariate values at site j, point p, visit r: Year, Type, Day (Julian day) and State

# number of observations for each site
Tj<-array(NA,j)
for (i in 1:j){
Tj[i]<-length(count.data$Pair2[glm.data$Pair2==pair[i]])
}

site<-sort(unique(count.data$Site))
j<-length(site)
years<-sort(unique(count.data$Year))
States<-sort(unique(count.data$State))

# matrix that will hold the counts n.jpr for each site (1 row per site)
Y<-matrix(NA,j,max(Tj))

# matrix that will hold the values for Year (1 for 2006, 2 for 2007, 3 for 2008): factor covariate with 2 levels
Year<-matrix(NA,j,max(Tj))

# matrix that will hold a 0 or 1 depending on whether Type = CONTROL or TREAT, respectively
Type<-matrix(NA,j,max(Tj))

# matrix that will hold the values for Julian day
Day<-matrix(NA,j,max(Tj))

# matrix that will hold the values for State
State<-matrix(NA,j,max(Tj))

# filling in the above values from
for (i in 1:j){

x<-which(count.data$Site==site[i])
l<-length(count.data$Count[x])
y<-count.data$Count[x]
    for (k in 1:l){

      Y[i,k]<-y[k]
      Year[i,k]<-ifelse(count.data$Year[x[k]]==2006,1,{ifelse(count.data$Year[x[k]]==2007,2,3)})
      Type[i,k]<-ifelse(count.data$Type[x[k]]=="TREAT",1,0)
      Day[i,k]<-count.data$jd[x[k]]
      State[i,k]<-which(States==count.data$State[x[k]])
      }
     }
sum(Y[which(is.na(Y)==F)])
#[1] 2545

########################################  set initial values ##############

# global hazard-rate model with scale and shape for L_y(\bmath{\theta}) (eqn (3)
scale0<-130
shape0<-2.5

#  count model parameters: intercept and random effect standard deviation for L_n(\bmath{\beta}|\bmath{\theta}) (eqn (6))
int0<- -13
std.ran0<-1

# the random effect coefficients b_j
b0<-rnorm(j,0,std.ran0)


#########################################################################
# setting up the matrices that will contain the paramter values;

# number of iterations
nt<-100000

# holds the values for detection function parameters for each iteration
# 15 colums due to 15 parameters in full model: hazard-rate det fct with covariates: year, type and state
det.param<-matrix(NA,nt+1,15)

# for an intercept only model:
det.param[1,]<-c(scale0, shape0, rep(0,13))

# holds the model id number for detection function for each iteration, 
det.model<-matrix(NA,nt+1,1)            # refers to det.list

# the matrix that will keep the parameter values for the density model
count.param<-matrix(NA,nt+1,(j+18))     

# filling in the initial values
count.param[1,1:(j+18)]<-c(int0,rep(0,16),std.ran0,b0)

# holds the model id number for density model for each iteration
count.model<-matrix(NA,nt+1,1)          # refers to count.list

############# proposal distributions
# proposal distributions for detection function parameters:

# 1. for main analysis and prior sensitivity analysis
det.prop.mean<-c(138.60, 3.00, 0.11, -0.15, 0.50, 0.42, 0.21, 0.70, 0.67, 0.64, 0.69, 0.61, 0.66, 0.03, 0.47)
det.prop.sd<-c(1.41, 0.84, rep(0.1,13))

# proposal distribution for the fixed effect density model parameters
# 1. for covey data
count.prop.mean<-c(-13.06, 0, 0.16, 0.08, 0.42, -0.01, 0, -0.71, -0.49, -1.16, -0.41, 0.01, -0.38, -1.36, 0.07, -1.05, 1.77)
count.prop.sd<-c(0.30, 0, 0.05, 0.05, 0.10,  0.01, 0, 0.24, 0.24, 0.23, 0.22, 0.20, 0.22, 0.23, 0.22, 0.23, 0.21)

################## picking the first model for detection function
det.model[1,1]<-5   # global hazard-rate: \bmath{\theta} = \{\sigma,\tau\}
cur.dmod<-det.model[1,1]
# holds the current det function parameters (vector \bmath{\theta}^t_m for model m)
rj.cursigs<-det.param[1,]


################## picking the first model for density model
# there are 16 models (all of them include random effect for Pair2):
# to pick the first model:
count.model[1]<-1          # intercept and random effect only model \bmath{\beta}=\{\beta_0,\sigma_b^2\}
cur.mod<-1
# holds the parameter values for the current density model (vector \bmath{\beta}^t for model m)
rj.curparam<-count.param[1,1:(j+18)]

# model identifier for detection function, model number in rows, parameters (y/n) in columns
det.list<-matrix(NA, 8, 15)
colnames(det.list)<-c("sig","sha","yea6","yea7","typ","sta","sta","sta","sta","sta","sta","sta","sta","sta","sta")
det.list[1,]<-c(1,1,1,1,0,rep(1,10))     # mcds with state and year
det.list[2,]<-c(1,1,0,0,1,rep(0,10))     # mcds with type
det.list[3,]<-c(1,1,0,0,0,rep(1,10))     # mcds with state
det.list[4,]<-c(1,1,0,0,1,rep(1,10))     # mcds with state and type
det.list[5,]<-c(1,1,rep(0,13))           # global hazard rate model
det.list[6,]<-c(1,1,rep(1,13))           # mcds with state year and type
det.list[7,]<-c(1,1,1,1,0,rep(0,10))     # mcds with year
det.list[8,]<-c(1,1,1,1,1,rep(0,10))     # mcds with year and type

# model identifier for density model, model number in rows, parameters (y/n) in columns
count.list<-matrix(0,16,17)
dimnames(count.list)[[2]]<-c("int","year6","year7","year8","typ","day","sta1","sta2","sta3","sta4","sta5","sta6","sta7","sta8","sta9","sta10","sta11")
count.list[1,c(1)]<-1
count.list[2,c(1,3:4)]<-1
count.list[3,c(1,5)]<-1
count.list[4,c(1,6)]<-1
count.list[5,c(1,3,4,5)]<-1
count.list[6,c(1,3,4,5,6)]<-1
count.list[7,c(1,3,4,6)]<-1
count.list[8,c(1,5,6)]<-1
count.list[9,c(1,8:17)]<-1
count.list[10,c(1,3,4,8:17)]<-1
count.list[11,c(1,5,8:17)]<-1
count.list[12,c(1,5,6,8:17)]<-1
count.list[13,c(1,3,4,5,8:17)]<-1
count.list[14,c(1,3,4,6,8:17)]<-1
count.list[15,c(1,3,4,5,6,8:17)]<-1
count.list[16,c(1,6,8:17)]<-1

################################# Likelihood functions ################################################

# hazard-rate function for point transects: \pi(y) * g(y|\bmath{\theta}) from eqn (2)
f.haz.function<-function(dis,sigma,shape) {
  f <- 2*pi*dis*(1-exp(-(dis/sigma)^(-shape)))
  return(f)
  }

############################### the priors ###########################################################

### detection function parameters
# for scale intercept
l.prior.sig<-function(sigm){
log.u.sig<-array(NA,length(sigm))
for (k in 1:length(sigm)){
log.u.sig[k]<-log(dunif(sigm[k],1,100000))                                
 }
return(sum(log.u.sig))}

# for shape
l.prior.sha<-function(shap){
log.u.sha<-array(NA,length(shap))
for (k in 1:length(shap)){
log.u<-log(dunif(shap[k],1,20))
ifelse(abs(log.u)==Inf,log.u.sha[k]<--100000,log.u.sha[k]<-log.u)}
return(sum(log.u.sha))}

# for scale coefficients
# for type covariate
l.prior.coeftyp<-function(coefsig){
lcs<-length(coefsig)
log.u.coefsig<-array(NA,lcs)
for (k in 1:lcs){
log.u<-log(dunif(coefsig[k],-3,3))          
ifelse(abs(log.u)==Inf,log.u.coefsig[k]<- -100000,log.u.coefsig[k]<-log.u)}
return(sum(log.u.coefsig))}

l.prior.coef<-function(coefsig){
lcs<-length(coefsig)
log.u.coefsig<-array(NA,lcs)
for (k in 1:lcs){
log.u<-log(dunif(coefsig[k],-2.5,2.5))
ifelse(abs(log.u)==Inf,log.u.coefsig[k]<- -100000,log.u.coefsig[k]<-log.u)}
return(sum(log.u.coefsig))}

### the priors for density model parameters
# intercept
l.prior.int<-function(int){
l.u.int<-log(dunif(int,-20,-7))
return(l.u.int)}

# year
l.prior.year<-function(yea){
l.u.yea1<-log(dunif(yea[1],-1,1));
l.u.yea2<-log(dunif(yea[2],-1,1));
l.u.yea<-sum(l.u.yea1,l.u.yea2)
if(abs(l.u.yea)==Inf){l.u.yea<- -100000}
return(l.u.yea)}

# 1 year
l.prior.year1<-function(yea){
l.u.yea<-log(dunif(yea[1],-1,1));
if(abs(l.u.yea)==Inf){l.u.yea<- -100000}
return(l.u.yea)}

# type
l.prior.type<-function(typ){
l.u.typ<-log(dunif(typ,0,3))
if(abs(l.u.typ)==Inf){l.u.typ<- -100000}
return(l.u.typ)}

# day
l.prior.day<-function(day){
l.u.day<-log(dunif(day,-0.1,0.1))
return(l.u.day)}

# state
l.prior.st<-function(st){
# last state is always zero (absorbed in intercept) so omitted
l.u.st1<-log(dunif(st[1],-3,3))
l.u.st2<-log(dunif(st[2],-3,3))
l.u.st3<-log(dunif(st[3],-3,3))
l.u.st4<-log(dunif(st[4],-3,3))
l.u.st5<-log(dunif(st[5],-3,3))
l.u.st6<-log(dunif(st[6],-3,3))
l.u.st7<-log(dunif(st[7],-3,3))
l.u.st8<-log(dunif(st[8],-3,3))
l.u.st9<-log(dunif(st[9],-3,3))
l.u.st10<-log(dunif(st[10],-3,3))
log.u.st<-sum(l.u.st1,l.u.st2,l.u.st3,l.u.st4,l.u.st5,l.u.st6,l.u.st7,l.u.st8,l.u.st9,l.u.st10)
if(abs(log.u.st)==Inf){log.u.st<- -100000 }
return(log.u.st)}

#  1 state
l.prior.1st<-function(st){
log.u.st<-log(dunif(st,-3,3))
if(abs(log.u.st)==Inf){log.u.st<- -100000}
return(log.u.st)}

# random effect standard deviation (std.ran)
l.prior.std.ran<-function(std.ran){
l.u.std.ran<-log(dunif(std.ran,0,2))
return(l.u.std.ran)}


########################## the posterior conditional distribution functions ######################

### integrated likelihood L_{n,y}(\bmath{\beta},\bmath{\theta}} from eqn (1)
# p combines all parameters \bmath{\theta}\bmath{\beta} and bj from eqns (1) - (7)
log.lik.fct<-function(p){
sig1<-p[1]               # det fct: scale intercept
sha2<-p[2]               # det fct: shape
sig.y<-c(p[3:4],0)       # det fct: year 2006,2007,2008 coef
sig.t<-p[5]              # det fct: type level CONTROL coef
sig.st<-c(p[6:15],0)     # det fct: states coef
int<-p[16]               # density intercept
yea<-p[17:19]            # density year 2006,2007,2008
typ<-p[20]               # density type level TREAT
day<-p[21]               # density day coef
st<-p[22:32]             # density state coef
std.ran<-p[33]           # density random effect standard deviation
b<-p[34:(j+33)]          # random effect coefficients

# 1. calculate the different scale parameters as function of parameters
sig.msyt<-matrix(NA,11,6)      # 11 states (rows), 3 years * 2 type levels (CONTROL,TREAT) (columns)
efa<-matrix(NA,11,6)

  for (strat in 1:11){        
   for (ty in 4:6){
   sig.msyt[strat,ty]<-sig1*exp(sig.st[strat]+ sig.y[ty-3])}}

   for (strat in 1:11){
    for (ty in 1:3){
    sig.msyt[strat,ty]<-sig.msyt[strat,ty+3]*exp(sig.t)}}

# 2. calculate the different effective areas as a function of covariates (using the scales from sig.msyt)
   for (strat in 1:11){
    for (ty in 1:6){
    efa[strat,ty]<-integrate(f.haz.function,0,500,sig.msyt[strat,ty],sha2)$value}}

# 3. calculate the f_e for each detection for det model likelihood component L_y(\bmath{\theta}) (eqn (3): exact distance data)
fe<-array(NA,length.d)
for (im in 1:length.d){
  ist<-which(states==covey.d$State[im])
  yst<-which(years==covey.d$Year[im])
  ifelse(covey.d$Type[im]=="CONTROL",multi<-0,multi<-3)
  norm.const<-efa[ist,multi+yst]    # normalising constant from denominator in eqn (2) equals the effective area for a given scale and shape parameter
  fe[im]<-log(f.haz.function(covey.d$Distance[im],sig.msyt[ist,multi+yst],sha2)/norm.const)
    }

# 4. arrange the offset from the current det model in a matrix for each observed count n.jpr in matrix Y (the \nu_{jpr} from eqn(6)):
offset<-matrix(NA,j,max(Tj))
for (is in 1:j){
 for (ts in 1:Tj[is]){
ifelse(Type[is,ts]==0,multi<-0,multi<-3)
offset[is,ts]<-efa[State[is,1],Year[is,ts]+multi]
}}

# 5. model L_n(\bmath{\beta}|\bmath{\theta}) from eqn (7)
 l.pois.y<-matrix(NA,j,max(Tj))  # matrix that will hold the Poisson likelihood for each observation n_jpr
 l.b.norm<-array(NA,j)           # vector that will hold the normal density for each random effect coefficients b_j
 lambda<-matrix(NA,j,max(Tj))    # matrix for storing the lambda_jpr from eqn (6)
   # for each site 
   for (ik in 1:j){
        # for each observation at that site
        for (k in 1:Tj[ik]){
        lambda[ik,k]<-exp(int + yea[Year[ik,k]] + Type[ik,k]*typ + day*Day[ik,k] + st[State[ik,k]] + b[ik] + log(offset[ik,k]))
        l.pois.y[ik,k]<-log(dpois(Y[ik,k],lambda[ik,k]))    # Poisson log-likelihood for each observation n_jpr in Y
        }
   l.b.norm[ik]<-log(dnorm(b[ik],0,std.ran))                # log of normal density for b_j
   }

   post<-sum(fe) + sum(l.pois.y[which(is.na(l.pois.y)==F)]) + sum(l.b.norm)
   return(post)
 }


################################# other function you will need
# this function matches a string of numbers against rows a matrix
# and will return the row number of matrix that matches the string exactly
match.function<-function(xmod,model.matrix){
is.match.or.not<-array(NA,length(model.matrix[,1]))
for (i in 1:length(model.matrix[,1])){
is.match.or.not[i]<-sum(xmod==model.matrix[i,])}
result<-which(is.match.or.not==length(model.matrix[1,]))
return(result)
}

i=2

#################################  the RJMCMC algorithm ######################################
# nt is the number of iterations and is set above
# row 1 is filled in with initial values for parameters and models
for (i in 2:nt){
print(i)

##################### RJ step : sequential proposals to add or delete covariates depending on whether they are in the model or not #####
# all models are considered equally likely, i.e. P(m|m') = P(m'|m) for all m' and m

############## the detection function  ########################

      # the current model
      cur.dmod<-det.model[i-1]
      curpa<-det.list[cur.dmod,]
      # setting the parameters for the new model equal to the current model
      newpa<-curpa
      rj.newsigs<-rj.cursigs

      # going through the list of coefficients to check whether to add or remove one
      # year coefficient
      ifelse(curpa[3]==0,{        # if year is not currently in the model, propose to add it
      newpa[3:4]<-1
      rj.newsigs[3:4]<-rnorm(2,det.prop.mean[3:4],det.prop.sd[3:4])   # draw random samples from proposal distributions
      # the numerator of eqn (11)
      num<-log.lik.fct(c(rj.newsigs,rj.curparam))+ l.prior.coef(rj.newsigs[3:4])
      # the denominator of eqn (11)
      den<-log.lik.fct(c(rj.cursigs,rj.curparam))+sum(log(dnorm(rj.newsigs[3:4],msyt.prop.mean[3:4],msyt.prop.sd[3:4])))
        } ,
        {
      newpa[3:4]<-0               # if year is in the current model, propose to delete it
      rj.newsigs[3:4]<-0
      # the numerator of eqn (11)
      num<-log.lik.fct(c(rj.newsigs,rj.curparam)) + sum(log(dnorm(rj.cursigs[3:4],msyt.prop.mean[3:4],msyt.prop.sd[3:4])))
      # the denominator of eqn (11)
      den<-log.lik.fct(c(rj.cursigs,rj.curparam))+ l.prior.coef(rj.cursigs[3:4])
        })
                   A<-min(1,exp(num-den))                   # proposed move is accepted with probability A
                   V<-runif(1)
                    ifelse(V<=A,{                           # if move is accepted change current values to new values
                       rj.cursigs<-rj.newsigs               
                       curpa<-newpa},
                                {                             
                       rj.newsigs<-rj.cursigs               # if move is rejected, reset everything to current
                       newpa<-curpa
                                })

      # type coefficient
      ifelse(curpa[5]==0,{        # if type is not in the current model, propose to add it
      newpa[5]<-1
      rj.newsigs[5]<-rnorm(1,det.prop.mean[5],det.prop.sd[5])  # draw a random sample from proposal distribution
      num<-log.lik.fct(c(rj.newsigs,rj.curparam))+ l.prior.coeftyp(rj.newsigs[5])
      den<-log.lik.fct(c(rj.cursigs,rj.curparam))+ sum(log(dnorm(rj.newsigs[5],msyt.prop.mean[5],msyt.prop.sd[5])))
        } ,
        {
      newpa[5]<-0                 # if type is in the current model, propose to delete it
      rj.newsigs[5]<-0
      num<-log.lik.fct(c(rj.newsigs,rj.curparam))+ sum(log(dnorm(rj.cursigs[5],msyt.prop.mean[5],msyt.prop.sd[5])))
      den<-log.lik.fct(c(rj.cursigs,rj.curparam))+ l.prior.coeftyp(rj.cursigs[5])
        })
                   A<-min(1,exp(num-den))
                   V<-runif(1)
                   ifelse(V<=A,{                             # if new model is accepted
                       rj.cursigs<-rj.newsigs                
                       curpa<-newpa},
                                {                            # if model is rejected, reset everything
                       rj.newsigs<-rj.cursigs
                       newpa<-curpa
                                })


      # for the state coefficients
      ifelse(curpa[6]==0,{
      newpa[6:15]<-1
      newpara<-which(newpa==1)
      rj.newsigs[6:15]<-rnorm(10,det.prop.mean[6:15],det.prop.sd[6:15])
      num<-log.lik.fct(c(rj.newsigs,rj.curparam)) + l.prior.coef(rj.newsigs[6:15])
      den<-log.lik.fct(c(rj.cursigs,rj.curparam)) + sum(log(dnorm(rj.newsigs[6:15],msyt.prop.mean[6:15],msyt.prop.sd[6:15])))
        } ,
        {
      newpa[6:15]<-0
      rj.newsigs[6:15]<-0
      num<-log.lik.fct(c(rj.newsigs,rj.curparam))+ sum(log(dnorm(rj.cursigs[6:15],msyt.prop.mean[6:15],msyt.prop.sd[6:15])))
      den<-log.lik.fct(c(rj.cursigs,rj.curparam))+ l.prior.coef(rj.cursigs[6:15])
        })
        den<-den.l+den.pro
        num<-num.l+num.pro

                   A<-min(1,exp(num-den))
                   V<-runif(1)
                    ifelse(V<=A,{                             
                       rj.cursigs<-rj.newsigs               
                       curpara<-newpara
                       curpa<-newpa},
                                {                             
                       rj.newsigs<-rj.cursigs
                       newpara<-curpara
                       newpa<-curpa
                                })

      # which model did we end up with 
      cur.dmod<-match.function(curpa,det.list)
# record the model selection for the det fct in det.model for the i'th iteration
det.model[i]<-cur.dmod

#################### RJ step for density model ##################################
      rj.newparam<-rj.curparam
      # the current model:
      cur.mod<-count.model[i-1]
      # the parameters in the current model
      cur.par<-count.list[cur.mod,]
      new.par<-cur.par

      # for the year coefficient
      ifelse (cur.par[3]==0,                                    # if year is not in current model, propose to add it
                 {new.par[3:4]<-1
                  # obtain a parameter value from the proposal distribution
                  rj.newparam[3:4]<-rnorm(2,count.prop.mean[3:4],count.prop.sd[3:4])
                    # test which model is better
                    # the numerator of eqn (11)
                    num<-log.lik.fct(c(rj.cursigs,rj.newparam)) + l.prior.year(rj.newparam[3:4]) 
                    # the denominator of eqn (11)
                    den<-log.lik.fct(c(rj.cursigs,rj.curparam)) + sum(log(dnorm(rj.newparam[3:4],count.prop.mean[3:4],count.prop.sd[3:4])))
                    A<-min(1,exp(num-den))
                    V<-runif(1)
                    ifelse(V<=A,{                             
                       rj.curparam<-rj.newparam               
                       cur.par<-new.par},
                                {                             
                       rj.newparam<-rj.curparam
                       new.par<-cur.par
                                })
                 } ,
                 {
                  rj.newparam[3:4]<-0                            # if year is in the current model, propose to delete it
                  new.par[3:4]<-0
                    num<-log.lik.fct(c(rj.cursigs,rj.newparam))  + sum(log(dnorm(rj.curparam[3:4],count.prop.mean[3:4],count.prop.sd[3:4])))
                    den<-log.lik.fct(c(rj.cursigs,rj.curparam))  + l.prior.year(rj.curparam[3:4]) #+ sum(log(dnorm(rj.newparam[new.p],count.prop.mean[new.p],count.prop.sd[new.p])))
                    A<-min(1,exp(num-den))
                    V<-runif(1)
                    ifelse(V<=A,{                             
                       rj.curparam<-rj.newparam              
                       cur.par<-new.par},
                                {                             
                       rj.newparam<-rj.curparam
                       new.par<-cur.par
                                })
                 } )
  # the type coefficient
              ifelse(cur.par[5]==0,                              # if type is not in current model, propose to add it
                 {new.par[5]<-1
                  rj.newparam[5]<-rnorm(1,count.prop.mean[5],count.prop.sd[5])
                    num<-log.lik.fct(c(rj.cursigs,rj.newparam)) + l.prior.type(rj.newparam[5]) 
                    den<-log.lik.fct(c(rj.cursigs,rj.curparam))  + sum(log(dnorm(rj.newparam[5],count.prop.mean[5],count.prop.sd[5])))
                    A<-min(1,exp(num-den))
                    V<-runif(1)
                    ifelse(V<=A,{                             
                       rj.curparam<-rj.newparam               
                       cur.par<-new.par},
                                {                             
                       rj.newparam<-rj.curparam
                       new.par<-cur.par
                                })
                 } ,
                 {
                  rj.newparam[5]<-0                              # if type is in the current model, propose to delete it
                  new.par[5]<-0
                    num<-log.lik.fct(c(rj.cursigs,rj.newparam))  + sum(log(dnorm(rj.curparam[5],count.prop.mean[5],count.prop.sd[5])))
                    den<-log.lik.fct(c(rj.cursigs,rj.curparam)) + l.prior.type(rj.curparam[5]) 
                    A<-min(1,exp(num-den))
                    V<-runif(1)
                    ifelse(V<=A,{                             
                       rj.curparam<-rj.newparam               
                       cur.par<-new.par},
                                {                             
                       rj.newparam<-rj.curparam
                       new.par<-cur.par
                                })
                 }  )
  # the day coefficient
              ifelse(cur.par[6]==0,                              # if Julian day is not in the current model, propose to add it
                 {new.par[6]<-1
                  rj.newparam[6]<-rnorm(1,count.prop.mean[6],count.prop.sd[6])
                    num<-log.lik.fct(c(rj.cursigs,rj.newparam)) + l.prior.day(rj.newparam[6]) 
                    den<-log.lik.fct(c(rj.cursigs,rj.curparam)) + sum(log(dnorm(rj.newparam[6],count.prop.mean[6],count.prop.sd[6])))
                    A<-min(1,exp(num-den))
                    V<-runif(1)
                    ifelse(V<=A,{                            
                       rj.curparam<-rj.newparam              
                       cur.par<-new.par},
                                {                            
                       rj.newparam<-rj.curparam
                       new.par<-cur.par
                                })
                 }  ,
                 {
                  rj.newparam[6]<-0                              # if Julian day is in the current model, propose to delete it
                  new.par[6]<-0
                    num<-log.lik.fct(c(rj.cursigs,rj.newparam)) + sum(log(dnorm(rj.curparam[6],count.prop.mean[6],count.prop.sd[6])))
                    den<-log.lik.fct(c(rj.cursigs,rj.curparam)) + l.prior.day(rj.curparam[6]) 
                    A<-min(1,exp(num-den))
                    V<-runif(1)
                    ifelse(V<=A,{                             
                       rj.curparam<-rj.newparam               
                       cur.par<-new.par},
                                {                             
                       rj.newparam<-rj.curparam
                       new.par<-cur.par
                                })
                 } )
  # the state coefficient
              ifelse(cur.par[8]==0,                             # if state is not in the current model, propose to add it
                 {new.par[8:17]<-1
                  for (f in 8:17){
                  rj.newparam[f]<-rnorm(1,count.prop.mean[f],count.prop.sd[f])}
                    num<-log.lik.fct(c(rj.cursigs,rj.newparam)) + l.prior.st(rj.newparam[8:17]) 
                    den<-log.lik.fct(c(rj.cursigs,rj.curparam)) + sum(log(dnorm(rj.newparam[8:17],count.prop.mean[8:17],count.prop.sd[8:17])))
                    A<-min(1,exp(rjc.num-rjc.den))
                    V<-runif(1)
                    ifelse(V<=A,{                             
                       rj.curparam<-rj.newparam               
                       cur.par<-new.par},
                                {                             
                       rj.newparam<-rj.curparam
                       new.par<-cur.par
                                })
                 } ,
                 {
                  rj.newparam[8:17]<-0                          # if state is in the current model, propose to delete it
                  new.par[8:17]<-0
                    num<-log.lik.fct(c(rj.cursigs,rj.newparam)) + sum(log(dnorm(rj.curparam[8:17],count.prop.mean[8:17],count.prop.sd[8:17])))
                    den<-log.lik.fct(c(rj.cursigs,rj.curparam)) + l.prior.st(rj.curparam[8:17]) 
                    A<-min(1,exp(rjc.num-rjc.den))
                    V<-runif(1)
                    ifelse(V<=A,{                             
                       rj.curparam<-rj.newparam               
                       cur.par<-new.par},
                                {                             
                       rj.newparam<-rj.curparam
                       new.par<-cur.par
                                })
                 }  )
# which model did we end up with:
cur.mod<-match.function(cur.par,count.list)
count.model[i]<-cur.mod

########################## Metropolis Hastings update ########################################################

########## updating the detection function parameters
mh.newsigs<-rj.cursigs
mh.cursigs<-rj.cursigs

              # for scale intercept
              u<-rnorm(1,0,3.5)                
              if((mh.cursigs[1]+u)>1){         # prevents scale intercept to become < 0
              mh.newsigs[1]<-mh.cursigs[1]+u
              # the numerator of eqn (8)
              num<-log.lik.fct(c(mh.newsigs,rj.curparam)) + l.prior.sig(mh.newsigs[1])
              # the denominator of eqn (8)
              den<-log.lik.fct(c(mh.cursigs,rj.curparam)) + l.prior.sig(mh.cursigs[1])
    A<-min(1,exp(num-den))
    V<-runif(1)
    ifelse(V<=A,mh.cursigs<-mh.newsigs,mh.newsigs<-mh.cursigs)    }
              # for shape                  
              u<-rnorm(1,0,0.2)
              mh.newsigs[2]<-mh.cursigs[2]+u
              num<-log.lik.fct(c(mh.newsigs,rj.curparam)) + l.prior.sha(mh.newsigs[2])
              den<-log.lik.fct(c(mh.cursigs,rj.curparam)) + l.prior.sha(mh.cursigs[2])
    A<-min(1,exp(num-den))
    V<-runif(1)
    ifelse(V<=A,mh.cursigs<-mh.newsigs,mh.newsigs<-mh.cursigs)

              # for the year coefficients in the scale parameter: levels 2006 and 2007
              if(mh.cursigs[3]!=0){
              u<-rnorm(1,0,0.12)
              mh.newsigs[3]<-mh.cursigs[3]+u
              num<-log.lik.fct(c(mh.newsigs,rj.curparam)) + l.prior.coef(mh.newsigs[3])
              den<-log.lik.fct(c(mh.cursigs,rj.curparam)) + l.prior.coef(mh.cursigs[3])
    A<-min(1,exp(num-den))
    V<-runif(1)
    ifelse(V<=A,mh.cursigs<-mh.newsigs,mh.newsigs<-mh.cursigs)

              u<-rnorm(1,0,0.12)
              mh.newsigs[4]<-mh.cursigs[4]+u
              num<-log.lik.fct(c(mh.newsigs,rj.curparam)) + l.prior.coef(mh.newsigs[4])
              den<-log.lik.fct(c(mh.cursigs,rj.curparam)) + l.prior.coef(mh.cursigs[4])
    A<-min(1,exp(num-den))
    V<-runif(1)
    ifelse(V<=A,mh.cursigs<-mh.newsigs,mh.newsigs<-mh.cursigs)
                                 }
              # for the type coefficient in the scale parameter: level CONTROL
              if(mh.cursigs[5]!=0){
              u<-rnorm(1,0,0.12)
              mh.newsigs[5]<-mh.cursigs[5]+u
              num<-log.lik.fct(c(mh.newsigs,rj.curparam)) + l.prior.coeftyp(mh.newsigs[5])
              den<-log.lik.fct(c(mh.cursigs,rj.curparam)) + l.prior.coeftyp(mh.cursigs[5])
    A<-min(1,exp(num-den))
    V<-runif(1)
    ifelse(V<=A,mh.cursigs<-mh.newsigs,mh.newsigs<-mh.cursigs)
                                 }

              # for the state coefficients in the scale parameter:
              if(mh.cursigs[6]!=0){
              for (ip in 6:15){
              #u<-rnorm(1,0,0.01)          # acc prob = 80-95%
              u<-rnorm(1,0,0.12)
              mh.newsigs[ip]<-mh.cursigs[ip]+u
              num<-log.lik.fct(c(mh.newsigs,rj.curparam)) + l.prior.coef(mh.newsigs[ip])
              den<-log.lik.fct(c(mh.cursigs,rj.curparam)) + l.prior.coef(mh.cursigs[ip])
    A<-min(1,exp(num-den))
    V<-runif(1)
    ifelse(V<=A,mh.cursigs<-mh.newsigs,mh.newsigs<-mh.cursigs)
    }
                                 }
    # fill in the new parameter values
    det.param[i,]<-mh.cursigs
    rj.cursigs<-mh.cursigs

######### updating the density model parameters
    curparam<-rj.curparam
    newparam<-rj.curparam

    # the intercept
               {u<-rnorm(1,0,0.08)                        
               newparam[1]<-curparam[1]+u
               num<-log.lik.fct(newparam) + l.prior.int(newparam[1])
               den<-log.lik.fct(curparam) + l.prior.int(curparam[1])
               A<-min(1,exp(num-den))
               V<-runif(1)
               ifelse(V<=A,curparam[1]<-newparam[1],newparam[1]<-curparam[1])
               }
    # the year coefficients:           
              if(curparam[3]!=0){
                       for (m in 3:4){
                       u<-rnorm(1,0,0.1)
                       newparam[m]<-curparam[m]+u
                       num<-log.lik.fct(c(rj.cursigs,newparam)) + l.prior.year1(newparam[m])
                       den<-log.lik.fct(c(rj.cursigs,curparam)) + l.prior.year1(curparam[m])
              A<-min(1,exp(num-den))
              V<-runif(1)
              ifelse(V<=A,curparam[m]<-newparam[m],newparam[m]<-curparam[m])
                       }}
    # the type coefficient
              if(curparam[5]!=0){
                       u<-rnorm(1,0,0.06)
                       newparam[5]<-curparam[5]+u
                       num<-log.lik.fct(c(rj.cursigs,newparam)) + l.prior.type(newparam[5])
                       den<-log.lik.fct(c(rj.cursigs,curparam)) + l.prior.type(curparam[5])
              A<-min(1,exp(num-den))
              V<-runif(1)
              ifelse(V<=A,curparam[5]<-newparam[5],newparam[5]<-curparam[5])
                       }
    # day coefficient
              if(curparam[6]!=0){
                       u<-rnorm(1,0,0.02)
                       newparam[6]<-curparam[6]+u
                       num<-log.lik.fct(c(rj.cursigs,newparam))    + l.prior.day(newparam[6])
                       den<-log.lik.fct(c(rj.cursigs,curparam))    + l.prior.day(curparam[6])
              A<-min(1,exp(num-den))
              V<-runif(1)
              ifelse(V<=A,curparam[6]<-newparam[6],newparam[6]<-curparam[6])
                       }
    # the state coefficients 
              ifelse(curparam[8]==0){
              for (m in 8:17){
                       u<-rnorm(1,0,0.25)
                       newparam[m]<-curparam[m]+u
                       num<-log.lik.fct(c(rj.cursigs,newparam))    + l.prior.1st(newparam[m])
                       den<-log.lik.fct(c(rj.cursigs,curparam))    + l.prior.1st(curparam[m])
              A<-min(1,exp(num-den))
              V<-runif(1)
              ifelse(V<=A,curparam[m]<-newparam[m],newparam[m]<-curparam[m])
              } }  
    # the random effect standard deviation
              {u<-max(rnorm(1,0,0.08),-newparam[18])    # cannot become 0 or less
              newparam[18]<-curparam[18]+u
              num<-log.ran.fct(c(rj.cursigs,newparam))    + l.prior.std.ran(newparam[18])
              den<-log.ran.fct(c(rj.cursigs,curparam))    + l.prior.std.ran(curparam[18])
              A<-min(1,exp(num-den))
              V<-runif(1)
              ifelse(V<=A,curparam[18]<-newparam[18],newparam[18]<-curparam[18])
              }
    # the random effects coefficients
    for (m in 19:(j+18)){
              u<-rnorm(1,0,0.4)
              newparam[m]<-curparam[m]+u
              num<-log.lik.fct(c(rj.cursigs,newparam))
              den<-log.lik.fct(c(rj.cursigs,curparam))
              A<-min(1,exp(num-den))
              V<-runif(1)
              ifelse(V<=A,curparam[m]<-newparam[m],newparam[m]<-curparam[m])
              }

    # saving the new parameter values of the density model in count.param
    count.param[i,]<-curparam[1:(j+18)]

    rj.curparam<-curparam
    rj.newparam<-curparam
  # saving the parameter matrices ever 1000 iterations
  if(!is.na(match(i,seq(0,200000,1000))==T)){
  save(det.model,file='det.model.RData')
  save(count.model,file='count.model.RData')
  save(det.param,file='msyt.param.RData')
  save(count.param,file='count.param.RData')
  }
  
  } # end of iteration



###### Alterations to this algorithm: 
# dis = distance (perpendicular for lines and radial for points, sigma = scale parameters, shape = shape parameter

# using a half-normal detection function for point transects for f(y); \pi(y) * g(y|\bmath{\theta}) from eqn (2) is given by
f.hn.function<-function(dis,sigma){
  f <- 2*pi*dis*exp(-dis^2/(2*sigma^2))
  return(f)
  }

# using a half-normal detection function for line transects for f*y); g(y|\bmath{\theta}) from eqn (2) is given by (\pi(y) can be ommitted for line transects)
f.hn.function <- function(dis,sigma) {
  f <- exp(-dis^2/(2*sigma^2))
  f
}

# using a hazard-rate function for point transects for f(y); \pi(y) * g(y|\bmath{\theta}) from eqn (2) is given by
f.haz.function<-function(dis,sigma,shape) {
  f <- 2*pi*dis*(1-exp(-(dis/sigma)^(-shape)))
  return(f)
  }

# using a hazard-rate function for line transects for f(y); g(y|\bmath{\theta}) from eqn (2) is given by (\pi(y) can be ommitted for line transects)
f.haz.function<-function(dis,sigma,shape) {
  f <- 1-exp(-(dis/sigma)^(-shape))
  return(f)
  }

# interval distance data: 
# the fi's for interval data using half-normal detection function for lines or points
# requires defining cutpoints for the intervals and the truncation distance w
f.haz.int<-function(cutpoint1,cutpoint2,sigma,shape){
 norm.const<-integrate(f.hn.function,0,w,sigma)$value
 int.prob<-integrate(f.haz.function,cutpoint1,cutpoint2,sigma,shape)$value/norm.const
 return(int.prob)
 }

# the fi's for interval data using hazard-rate detection function for lines or points
# requires defining cutpoints for the intervals and the truncation distance w
f.haz.int<-function(cutpoint1,cutpoint2,sigma,shape){
 norm.const<-integrate(f.haz.function,0,w,sigma,shape)$value
 int.prob<-integrate(f.haz.function,cutpoint1,cutpoint2,sigma,shape)$value/norm.const
 return(int.prob)
 }
