#This code generates a table as given for the M-bridge data study. This code uses a dummy real life data 
#emulating the M-bridge data, since the original data cannot be shared due to privacy. The result is displayed in
#the console after execution of the code.

#loading library
library(parallel)
library(dplyr)
library(ggplot2)

#Reading the real life example data set from the working directory
data = read.csv("real_life_example_data.csv")
data = data[,-1]
data = as.data.frame(data)

resp_est_prob = function(data,t1,t2){
  pattr = which(data[,1] == t1 & data[,2] == t2) #checking number of patients in intervention
  detr = which(data[,1] == t1)
  if(length(pattr) == 0){
    p_tr = 0}else{
      p_tr = length(pattr)/length(detr)
    }
}

wash_out_suc = function(n){
  s = sum(data$bingefin[n])
  return(s)
}

prob_est_func = function(n,s){
  if(length(n) == 0 | s == 0){
    p = 0.01
  }else{
    p = s/length(n)
  }
  if(p == 1){
    p = 0.99
  }else{
    p = p
  }
  return(p)
}

vec_tau_a = NULL
vec_tau_ac = NULL
vec_tau_be = NULL

g_a_est = resp_est_prob(data,"A","A")
g_b_est = resp_est_prob(data,"B","B")

gp_a = g_a_est
gp_b = g_b_est

waiting_lot = 60

n_AA = which(data$First.Treatment == "A" & data$Second.Treatment == "A")[1:(round(waiting_lot*0.5)-round(waiting_lot*0.5*(1-g_a_est)*0.5)-round(waiting_lot*0.5*(1-g_a_est)*0.5))]
n_AC = which(data$First.Treatment == "A" & data$Second.Treatment == "C")[1:(round(waiting_lot*0.5*(1-g_a_est)*0.5))]
n_AD = which(data$First.Treatment == "A" & data$Second.Treatment == "D")[1:(round(waiting_lot*0.5*(1-g_a_est)*0.5))]
n_BB = which(data$First.Treatment == "B" & data$Second.Treatment == "B")[1:(round(waiting_lot*0.5)-round(waiting_lot*0.5*(1-g_b_est)*0.5)-round(waiting_lot*0.5*(1-g_b_est)*0.5))]
n_BE = which(data$First.Treatment == "B" & data$Second.Treatment == "E")[1:(round(waiting_lot*0.5*(1-g_b_est)*0.5))]
n_BF = which(data$First.Treatment == "B" & data$Second.Treatment == "F")[1:(round(waiting_lot*0.5*(1-g_b_est)*0.5))]

n_A = c(n_AA,n_AC,n_AD)
n_B = c(n_BB,n_BE,n_BF)

n = max(n_A,n_B)

s_AA = wash_out_suc(n_AA)
s_AC = wash_out_suc(n_AC)
s_AD = wash_out_suc(n_AD)
s_BB = wash_out_suc(n_BB)
s_BE = wash_out_suc(n_BE)
s_BF = wash_out_suc(n_BF)

p_aa = prob_est_func(n_AA,s_AA)
p_ac = prob_est_func(n_AC,s_AC)
p_ad = prob_est_func(n_AD,s_AD)
p_bb = prob_est_func(n_BB,s_BB)
p_be = prob_est_func(n_BE,s_BE)
p_bf = prob_est_func(n_BF,s_BF)

p_a = (g_a_est*p_aa)+((1-g_a_est)*((p_ac)^1.5)/((p_ac)^0.5+(p_ad)^0.5))+((1-g_a_est)*((p_ad)^1.5)/((p_ac)^0.5+(p_ad)^0.5))
p_b = (g_b_est*p_bb)+((1-g_b_est)*((p_be)^1.5)/((p_be)^0.5+(p_bf)^0.5))+((1-g_b_est)*((p_bf)^1.5)/((p_be)^0.5+(p_bf)^0.5))

#calculating the first stage optimal allocation probability
tau_a = (sqrt(p_a/p_b))
alloc_a = tau_a/(1+tau_a)
vec_tau_a = c(vec_tau_a,tau_a)
alloc_a = tau_a/(1+tau_a)

#calculating the second stage optimal allocation probability along A
tau_ac = (sqrt(p_ac/p_ad))
vec_tau_ac = c(vec_tau_ac,tau_ac)
alloc_ac = tau_ac/(1+tau_ac)

#calculating the second stage optimal allocation probability along B
tau_be = (sqrt(p_be/p_bf))
vec_tau_be = c(vec_tau_be,tau_be)
alloc_be = tau_be/(1+tau_be)

inter_data = NULL
inter_data = data[-c(n_A,n_B),]
inter_data = as.data.frame(inter_data)

to_prin_row = NULL
to_prin = NULL
to_prin_row = c(length(n_A)+length(n_B),length(n_A),length(n_B),length(n_AA),length(n_AC),length(n_AD),length(n_BB),length(n_BE),length(n_BF),gp_a,gp_b,p_aa,p_ac,p_ad,p_bb,p_be,p_bf,tau_a,tau_ac,tau_be,alloc_a,alloc_ac,alloc_be)
to_prin = rbind(to_prin,to_prin_row)
to_prin_row = NULL

while(1>0){
  pat_lot = 1
  lot_AA = NULL
  lot_AC = NULL
  lot_AD = NULL
  lot_BB = NULL
  lot_BE = NULL
  lot_BF = NULL
  
  ft_rand = rbinom(1,1,alloc_a)
  res_a = rbinom(1,1,g_a_est)
  sc_rand_a = rbinom(1,1,alloc_ac)
  res_b = rbinom(1,1,g_b_est)
  sc_rand_b = rbinom(1,1,alloc_be)
  
  if(length(which(inter_data$First.Treatment == "A" & inter_data$Second.Treatment == "A")) == 0)
    break
  if(length(which(inter_data$First.Treatment == "A" & inter_data$Second.Treatment == "C")) == 0)
    break
  if(length(which(inter_data$First.Treatment == "A" & inter_data$Second.Treatment == "D")) == 0)
    break
  if(length(which(inter_data$First.Treatment == "B" & inter_data$Second.Treatment == "B")) == 0)
    break
  if(length(which(inter_data$First.Treatment == "B" & inter_data$Second.Treatment == "E")) == 0)
    break
  if(length(which(inter_data$First.Treatment == "B" & inter_data$Second.Treatment == "F")) == 0)
    break
  
  
  if(ft_rand == 1 & res_a == 1){
    lot_AA = which(inter_data$First.Treatment == "A" & inter_data$Second.Treatment == "A")[1]
  }
  if(ft_rand == 1 & res_a == 0){
    if(sc_rand_a == 1){
      lot_AC = which(inter_data$First.Treatment == "A" & inter_data$Second.Treatment == "C")[1]
    }
    if(sc_rand_a == 0){
      lot_AD = which(inter_data$First.Treatment == "A" & inter_data$Second.Treatment == "D")[1]
    }
  }
  n_AA = c(n_AA,lot_AA)
  n_AC = c(n_AC,lot_AC)
  n_AD = c(n_AD,lot_AD)
  
  if(ft_rand == 0 & res_b == 1){
    lot_BB = which(inter_data$First.Treatment == "B" & inter_data$Second.Treatment == "B")[1]
  }
  if(ft_rand == 0 & res_b == 0){
    if(sc_rand_b == 1){
      lot_BE = which(inter_data$First.Treatment == "B" & inter_data$Second.Treatment == "E")[1]
    }
    if(sc_rand_b == 0){
      lot_BF = which(inter_data$First.Treatment == "B" & inter_data$Second.Treatment == "F")[1]
    }
  }
  
  n_BB = c(n_BB,lot_BB)
  n_BE = c(n_BE,lot_BE)
  n_BF = c(n_BF,lot_BF)
  
  n_A = c(n_AA,n_AC,n_AD)
  n_B = c(n_BB,n_BE,n_BF)
  
  drop_pat = c(lot_AA,lot_AC,lot_AD,lot_BB,lot_BE,lot_BF)
  
  s_AA = s_AA+sum(inter_data$Final.Outcome[lot_AA])
  s_AC = s_AC+sum(inter_data$Final.Outcome[lot_AC])
  s_AD = s_AD+sum(inter_data$Final.Outcome[lot_AD])
  s_BB = s_BB+sum(inter_data$Final.Outcome[lot_BB])
  s_BE = s_BE+sum(inter_data$Final.Outcome[lot_BE])
  s_BF = s_BF+sum(inter_data$Final.Outcome[lot_BF])
  
  p_aa = prob_est_func(n_AA,s_AA)
  p_ac = prob_est_func(n_AC,s_AC)
  p_ad = prob_est_func(n_AD,s_AD)
  p_bb = prob_est_func(n_BB,s_BB)
  p_be = prob_est_func(n_BE,s_BE)
  p_bf = prob_est_func(n_BF,s_BF)
  
  gp_a = length(n_AA)/length(n_A)
  gp_b = length(n_BB)/length(n_B)
  
  #calculating first stage estimated success probabilities
  p_a = (g_a_est*p_aa)+((1-g_a_est)*((p_ac)^1.5)/((p_ac)^0.5+(p_ad)^0.5))+((1-g_a_est)*((p_ad)^1.5)/((p_ac)^0.5+(p_ad)^0.5))
  p_b = (g_b_est*p_bb)+((1-g_b_est)*((p_be)^1.5)/((p_be)^0.5+(p_bf)^0.5))+((1-g_b_est)*((p_bf)^1.5)/((p_be)^0.5+(p_bf)^0.5))
  
  
  # #calculating the derivatives for treatment A
  # f_d_p_ac = (1-gp_a)*((3*p_ac*(p_ad)^0.5+2*(p_ac)^1.5-(p_ad)^1.5)/(2*(p_ac)^0.5*((p_ac)^0.5+(p_ad)^0.5)^2))
  # f_d_p_ad = (1-gp_a)*((3*p_ad*(p_ac)^0.5+2*(p_ad)^1.5-(p_ac)^1.5)/(2*(p_ad)^0.5*((p_ac)^0.5+(p_ad)^0.5)^2))
  # #f_d_q_ac = (1-g_a)*((-3*p_ac*(p_ad)^0.5-2*(p_ac)^1.5+(p_ad)^1.5)/(2*(p_ac)^0.5*((p_ac)^0.5+(p_ad)^0.5)^2))
  # #f_d_q_ac = (1-g_a)*((-3*p_ad*(p_ac)^0.5-2*(p_ad)^1.5+(p_ac)^1.5)/(2*(p_ad)^0.5*((p_ac)^0.5+(p_ad)^0.5)^2))
  # s_d_p_ac = (1-gp_a)*(((p_ac)^2*(p_ad)^0.5+6*p_ac*(p_ad)^1.5+4*(p_ad)^2*(p_ac)^0.5+4*p_ad*(p_ac)^1.5+(p_ad)^2.5)/(4*(p_ac)^1.5*((p_ac)^0.5+(p_ad)^0.5)^4))
  # s_d_p_ad = (1-gp_a)*(((p_ad)^2*(p_ac)^0.5+6*p_ad*(p_ac)^1.5+4*(p_ac)^2*(p_ad)^0.5+4*p_ac*(p_ad)^1.5+(p_ac)^2.5)/(4*(p_ad)^1.5*((p_ac)^0.5+(p_ad)^0.5)^4))
  # #s_d_q_ac = (g_a-1)*(((p_ac)^2*(p_ad)^0.5+6*p_ac*(p_ad)^1.5+4*(p_ad)^2*(p_ac)^0.5+4*p_ad*(p_ac)^1.5+(p_ad)^2.5)/(4*(p_ac)^1.5*((p_ac)^0.5+(p_ad)^0.5)^4))
  # #s_d_q_ad = (g_a-1)*(((p_ad)^2*(p_ac)^0.5+6*p_ad*(p_ac)^1.5+4*(p_ac)^2*(p_ad)^0.5+4*p_ac*(p_ad)^1.5+(p_ac)^2.5)/(4*(p_ad)^1.5*((p_ac)^0.5+(p_ad)^0.5)^4))
  # 
  # #calculating the derivatives for treatment B
  # f_d_p_be = (1-gp_b)*((3*p_be*(p_bf)^0.5+2*(p_be)^1.5-(p_bf)^1.5)/(2*(p_be)^0.5*((p_be)^0.5+(p_bf)^0.5)^2))
  # f_d_p_bf = (1-gp_b)*((3*p_bf*(p_be)^0.5+2*(p_bf)^1.5-(p_be)^1.5)/(2*(p_bf)^0.5*((p_be)^0.5+(p_bf)^0.5)^2))
  # #f_d_q_ac = (1-g_a)*((-3*p_ac*(p_ad)^0.5-2*(p_ac)^1.5+(p_ad)^1.5)/(2*(p_ac)^0.5*((p_ac)^0.5+(p_ad)^0.5)^2))
  # #f_d_q_ac = (1-g_a)*((-3*p_ad*(p_ac)^0.5-2*(p_ad)^1.5+(p_ac)^1.5)/(2*(p_ad)^0.5*((p_ac)^0.5+(p_ad)^0.5)^2))
  # s_d_p_be = (1-gp_b)*(((p_be)^2*(p_bf)^0.5+6*p_be*(p_bf)^1.5+4*(p_bf)^2*(p_be)^0.5+4*p_bf*(p_be)^1.5+(p_bf)^2.5)/(4*(p_be)^1.5*((p_be)^0.5+(p_bf)^0.5)^4))
  # s_d_p_bf = (1-gp_b)*(((p_bf)^2*(p_be)^0.5+6*p_bf*(p_be)^1.5+4*(p_be)^2*(p_bf)^0.5+4*p_be*(p_bf)^1.5+(p_be)^2.5)/(4*(p_bf)^1.5*((p_be)^0.5+(p_bf)^0.5)^4))
  # #s_d_q_ac = (g_a-1)*(((p_ac)^2*(p_ad)^0.5+6*p_ac*(p_ad)^1.5+4*(p_ad)^2*(p_ac)^0.5+4*p_ad*(p_ac)^1.5+(p_ad)^2.5)/(4*(p_ac)^1.5*((p_ac)^0.5+(p_ad)^0.5)^4))
  # #s_d_q_ad = (g_a-1)*(((p_ad)^2*(p_ac)^0.5+6*p_ad*(p_ac)^1.5+4*(p_ac)^2*(p_ad)^0.5+4*p_ac*(p_ad)^1.5+(p_ac)^2.5)/(4*(p_ad)^1.5*((p_ac)^0.5+(p_ad)^0.5)^4))
  
  #calculating the first stage optimal allocation probability
  tau_a = (sqrt(p_a/p_b))
  alloc_a = tau_a/(1+tau_a)
  vec_tau_a = c(vec_tau_a,tau_a)
  alloc_a = tau_a/(1+tau_a)
  
  #calculating the second stage optimal allocation probability along A
  tau_ac = (sqrt(p_ac/p_ad))
  vec_tau_ac = c(vec_tau_ac,tau_ac)
  alloc_ac = tau_ac/(1+tau_ac)
  
  #calculating the second stage optimal allocation probability along B
  tau_be = (sqrt(p_be/p_bf))
  vec_tau_be = c(vec_tau_be,tau_be)
  alloc_be = tau_be/(1+tau_be)
  
  inter_data = inter_data[-drop_pat,]
  
  #n_A = c(n_A,which(inter_data$Group.x == "Experimental Early" & inter_data$Heavy_Drinker.x == "")[1:(pat_lot*p_a)])
  #n_B = c(n_B,which(inter_data$Group.x == "Experimental Late" & inter_data$Heavy_Drinker.x == "")[1:(pat_lot*p_b)])
  
  to_prin_row = c(length(n_A)+length(n_B),length(n_A),length(n_B),length(n_AA),length(n_AC),length(n_AD),length(n_BB),length(n_BE),length(n_BF),gp_a,gp_b,p_aa,p_ac,p_ad,p_bb,p_be,p_bf,tau_a,tau_ac,tau_be,alloc_a,alloc_ac,alloc_be)
  
  
  # print(length(which(inter_data$Group.y == "Experimental Early" & inter_data$Heavy_Drinker.y == "")))
  # print(sum(inter_data$bingefin[which(inter_data$Group.y == "Experimental Early" & inter_data$Heavy_Drinker.y == "")]))
  # print(length(which(inter_data$Group.y == "Experimental Early" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Mbridge.y == "TRUE")))
  # print(sum(inter_data$bingefin[which(inter_data$Group.y == "Experimental Early" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Mbridge.y == "TRUE")]))
  # print(length(which(inter_data$Group.y == "Experimental Early" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Auto_Email.y == "TRUE")))
  # print(sum(inter_data$bingefin[which(inter_data$Group.y == "Experimental Early" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Auto_Email.y == "TRUE")]))
  # 
  # print(length(which(inter_data$Group.y == "Experimental Late" & inter_data$Heavy_Drinker.y == "")))
  # print(sum(inter_data$bingefin[which(inter_data$Group.y == "Experimental Late" & inter_data$Heavy_Drinker.y == "")]))
  # print(length(which(inter_data$Group.y == "Experimental Late" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Mbridge.y == "TRUE")))
  # print(sum(inter_data$bingefin[which(inter_data$Group.y == "Experimental Late" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Mbridge.y == "TRUE")]))
  # print(length(which(inter_data$Group.y == "Experimental Late" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Auto_Email.y == "TRUE")))
  # print(sum(inter_data$bingefin[which(inter_data$Group.y == "Experimental Late" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Auto_Email.y == "TRUE")]))
  # 
  # 
  # print(to_prin_row)
  to_prin = rbind(to_prin,to_prin_row)
  to_prin_row = NULL
}

#View(to_prin)



to_prin = as.data.frame(to_prin)
to_prin[,c(10:23)] = round(to_prin[,c(10:23)],2)

rest_AA = length(which(inter_data$First.Treatment == "A" & inter_data$Second.Treatment == "A"))
rest_AC = length(which(inter_data$First.Treatment == "A" & inter_data$Second.Treatment == "C"))
rest_AD = length(which(inter_data$First.Treatment == "A" & inter_data$Second.Treatment == "D"))
rest_BB = length(which(inter_data$First.Treatment == "B" & inter_data$Second.Treatment == "B"))
rest_BE = length(which(inter_data$First.Treatment == "B" & inter_data$Second.Treatment == "E"))
rest_BF = length(which(inter_data$First.Treatment == "B" & inter_data$Second.Treatment == "F"))
rest_total = length(which(inter_data$First.Treatment == "A"))+length(which(inter_data$First.Treatment == "B"))

if(rest_AA == 0){
  failure_AA = 0
}else{
  failure_AA = (rest_AA-(length(which(inter_data$Group.y == "Experimental Early" & inter_data$Heavy_Drinker.y == "" & inter_data$bingefin == 0))))/(rest_AA)
}
if(rest_AC == 0){
  failure_Ac = 0
}else{
  failure_AC = (rest_AC-(length(which(inter_data$Group.y == "Experimental Early" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Mbridge.y == "TRUE" & inter_data$bingefin == 0))))/(rest_AC)
}
if(rest_AD == 0){
  failure_AD = 0
}else{
  failure_AD = (rest_AD-(length(which(inter_data$Group.y == "Experimental Early" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Auto_Email.y == "TRUE" & inter_data$bingefin == 0))))/(rest_AD)
}
if(rest_BB == 0){
  failure_BB = 0
}else{
  failure_BB = (rest_BB-(length(which(inter_data$Group.y == "Experimental Late" & inter_data$Heavy_Drinker.y == "" & inter_data$bingefin == 0))))/(rest_BB)
}
if(rest_BE == 0){
  failure_BE = 0
}else{
  failure_BE = (rest_BE-(length(which(inter_data$Group.y == "Experimental Late" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Mbridge.y == "TRUE" & inter_data$bingefin == 0))))/(rest_BE)
}
if(rest_BF == 0){
  failure_BF = 0
}else{
  failure_BF = (rest_BF-(length(which(inter_data$Group.y == "Experimental Late" & inter_data$Heavy_Drinker.y == "Yes" & inter_data$Intervention_Auto_Email.y == "TRUE" & inter_data$bingefin == 0))))/(rest_BF)
}

failure_rest_d1 = g_a_est*failure_AA+(1-g_a_est)*failure_AC
failure_rest_d2 = g_a_est*failure_AA+(1-g_a_est)*failure_AD
failure_rest_d3 = g_b_est*failure_BB+(1-g_b_est)*failure_BE
failure_rest_d4 = g_b_est*failure_BB+(1-g_b_est)*failure_BF
failure_rest_total = (length(which(inter_data$Group.y == "Experimental Early" & inter_data$bingefin == 0))+length(which(inter_data$Group.y == "Experimental Late" & inter_data$bingefin == 0)))/(rest_total)


rbind(cbind(length(n_AA)+length(n_AC),g_a_est*((length(n_AA)-s_AA)/length(n_AA))+(1-g_a_est)*((length(n_AC)-s_AC)/length(n_AC)),(rest_AA+rest_AC),failure_rest_d1),
      cbind(length(n_AA)+length(n_AD),g_a_est*((length(n_AA)-s_AA)/length(n_AA))+(1-g_a_est)*((length(n_AD)-s_AD)/length(n_AD)),(rest_AA+rest_AD),failure_rest_d2),
      cbind(length(n_BB)+length(n_BE),g_b_est*((length(n_BB)-s_BB)/length(n_BB))+(1-g_b_est)*((length(n_BE)-s_BE)/length(n_BE)),(rest_BB+rest_BE),failure_rest_d3),
      cbind(length(n_BB)+length(n_BF),g_b_est*((length(n_BB)-s_BB)/length(n_BB))+(1-g_b_est)*((length(n_BF)-s_BF)/length(n_BF)),(rest_BB+rest_BF),failure_rest_d4),
      cbind(length(n_A)+length(n_B),(length(n_A)+length(n_B)-(s_AA+s_AC+s_AD+s_BB+s_BE+s_BF))/(length(n_A)+length(n_B)),rest_total,failure_rest_total))

prop_vec = c(g_a_est*(s_AA/length(n_AA))+(1-g_a_est)*(s_AC/length(n_AC)),g_a_est*(s_AA/length(n_AA))+(1-g_a_est)*(s_AD/length(n_AD)),g_b_est*(s_BB/length(n_BB))+(1-g_b_est)*(s_BE/length(n_BE)),g_b_est*(s_BB/length(n_BB))+(1-g_b_est)*(s_BF/length(n_BF)))
prop_vec_1 = c(g_a_est*((length(n_AA)-s_AA)/length(n_AA))+(1-g_a_est)*((length(n_AC)-s_AC)/length(n_AC)),g_a_est*((length(n_AA)-s_AA)/length(n_AA))+(1-g_a_est)*((length(n_AD)-s_AD)/length(n_AD)),g_b_est*((length(n_BB)-s_BB)/length(n_BB))+(1-g_b_est)*((length(n_BE)-s_BE)/length(n_BE)),g_b_est*((length(n_BB)-s_BB)/length(n_BB))+(1-g_b_est)*((length(n_BF)-s_BF)/length(n_BF)))
patient_vec = c(length(n_AA)+length(n_AC),length(n_AA)+length(n_AD),length(n_BB)+length(n_BE),length(n_BB)+length(n_BF))

final_vec = c(length(n_AA)+length(n_AC),length(n_AA)+length(n_AD),length(n_BB)+length(n_BE),length(n_BB)+length(n_BF),prop_vec_1)
final_vec = as.data.frame(t(final_vec))
names = c("Patients_d1","Patients_d2","Patients_d3","Patients_d4","failure_prop_d1","failure_prop_d2","failure_prop_d3","failure_prop_d4")
colnames(final_vec) = names

print(final_vec)
