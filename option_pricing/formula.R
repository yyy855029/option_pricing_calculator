library(tidyverse)
# grid.arrange function
library(gridExtra)


# tab1
# Black Scholes Model
bs_price<-function(S,K,r,maturity,sig){
  if(maturity>0){
    d1<-(log(S/K)+(r+0.5*sig^2)*maturity)/(sig*sqrt(maturity))
    d2<-d1-sig*sqrt(maturity)
    N_d1<-pnorm(d1)
    n_d1<-dnorm(d1)
    N_d2<-pnorm(d2)
    n_d2<-dnorm(d2)
  }
  else{
    if(S>=K){
      N_d1<-1
      n_d1<-0
    }
    else{
      N_d1<-0
      n_d1<-0
    }
    N_d2<-N_d1
    n_d2<-n_d1
  }
  
  # Price
  Call<-S*N_d1-K*exp(-r*maturity)*N_d2 
  Put<-K*exp(-r*maturity)*(1-N_d2)-S*(1-N_d1)
  
  return(list(Call,Put,N_d1,N_d2,n_d1,n_d2))
}

make_bs_price_df<-function(S,K,r,maturity,sig){
  price_list<-bs_price(S,K,r,maturity,sig)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  N_d1<-price_list[[3]]
  N_d2<-price_list[[4]]
  n_d1<-price_list[[5]]
  n_d2<-price_list[[6]]
  
  N_md1<-1-N_d1
  N_md2<-1-N_d2
  
  # Greek
  delta_C<-N_d1
  gamma_C<-n_d1/(S*sig*sqrt(maturity))
  vega_C<-S*n_d1*sqrt(T)
  rho_C<-K*maturity*exp(-r*maturity)*N_d2
  theta_C<--(S*n_d1*sig/(2*sqrt(T))+K*r*exp(-r*maturity)*N_d2)
  
  delta_P<--N_md1
  gamma_P<-gamma_C
  vega_P<-vega_C
  rho_P<--K*maturity*exp(-r*maturity)*N_md2
  theta_P<--(S*n_d1*sig/(2*sqrt(T))-K*r*exp(-r*maturity)*N_md2)
  
  # Create dataframe
  bs_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(bs_cp_df)<-c('Call Price','Put Price')
  bs_N_d_df<-data.frame('N(d1)'=N_d1,'N(d2)'=N_d2)
  colnames(bs_N_d_df)<-c('N (d1)','N (d2)')
  bs_greek_df<-data.frame(Delta=c(delta_C,delta_P),Gamma=c(gamma_C,gamma_P),Vega=c(vega_C,vega_P),
                          Rho=c(rho_C,rho_P),Theta=c(theta_C,theta_P))
  rownames(bs_greek_df)<-c('Call','Put')
  
  return(list(bs_cp_df,bs_N_d_df,bs_greek_df))
}


# Binomial Tree Model (European)
bs_bt_euro_price<-function(S,K,r,maturity,sig,tree_n){
  t<-maturity/tree_n
  u<-exp(sig*sqrt(t))
  d<-1/u
  p<-(exp(r*t)-d)/(u-d)
  q<-1-p
  FV_C<-0
  FV_P<-0
  
  for(i in seq(tree_n+1)){
    FV_C<-FV_C+choose(tree_n,i-1)*((p^(tree_n-i+1))*(q^(i-1)))*pmax(S*((u^(tree_n-i+1))*(d^(i-1)))-K,0)
    FV_P<-FV_P+choose(tree_n,i-1)*((p^(tree_n-i+1))*(q^(i-1)))*pmax(K-S*((u^(tree_n-i+1))*(d^(i-1))),0)
  }
  
  # Price
  Call<-exp(-r*maturity)*FV_C
  Put<-exp(-r*maturity)*FV_P
  
  return(list(Call,Put))
}

make_bs_bt_euro_price_df<-function(S,K,r,maturity,sig,tree_n){
  price_list<-bs_bt_euro_price(S,K,r,maturity,sig,tree_n)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  
  # Create dataframe
  bt_euro_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(bt_euro_cp_df)<-c('Call Price','Put Price')
  
  return(bt_euro_cp_df)
}


# Binomial Tree Model (American)
bs_bt_amer_price<-function(S,K,r,maturity,sig,tree_n){
  t<-maturity/tree_n
  u<-exp(sig*sqrt(t))
  d<-1/u
  p<-(exp(r*t)-d)/(u-d)
  q<-1-p
  discount<-exp(-r*t)
  
  lattice_C<-matrix(0L,nrow=tree_n+1,ncol=tree_n+1) 
  lattice_C[1,1]<-S
  
  # i:第幾層
  for(i in seq(tree_n)){
    # j:在i層中第幾個
    for(j in seq(i)){
      lattice_C[i+1,j+1]<-u*lattice_C[i,j]
    }
    lattice_C[i+1,1]<-d*lattice_C[i,1]
  }
  
  lattice_P<-lattice_C
  
  for(i in seq(tree_n+1,2)){
    for(j in seq(i,2)){
      if(i==tree_n+1){
        euro_C<-discount*(p*pmax(lattice_C[i,j]-K,0)+q*pmax(lattice_C[i,j-1]-K,0))
        euro_P<-discount*(p*pmax(K-lattice_P[i,j],0)+q*pmax(K-lattice_P[i,j-1],0))
      }
      else{
        euro_C<-discount*(p*lattice_C[i,j]+q*lattice_C[i,j-1])
        euro_P<-discount*(p*lattice_P[i,j]+q*lattice_P[i,j-1])
      }
      exer_C<-pmax(lattice_C[i-1,j-1]-K,0)
      exer_P<-pmax(K-lattice_P[i-1,j-1],0)
      lattice_C[i-1,j-1]<-pmax(euro_C,exer_C)
      lattice_P[i-1,j-1]<-pmax(euro_P,exer_P)
    }
  }
  
  # Price
  Call<-lattice_C[1,1]
  Put<-lattice_P[1,1]
  
  return(list(Call,Put,u,d,p))
  
}

make_bs_bt_amer_price_df<-function(S,K,r,maturity,sig,tree_n){
  price_list<-bs_bt_amer_price(S,K,r,maturity,sig,tree_n)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  u<-price_list[[3]]
  d<-price_list[[4]]
  p<-price_list[[5]]
  
  # Create dataframe
  bt_amer_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(bt_amer_cp_df)<-c('Call Price','Put Price')
  bt_arg_df<-data.frame('Up'=u,'Down'=d,'Up Probability'=p)
  colnames(bt_arg_df)<-c('Up','Down','Up Probability')
  
  return(list(bt_amer_cp_df,bt_arg_df))
}


# Monte Carlo Simulation
bs_mc_price<-function(S,K,r,maturity,sig,sim_n){
  norm_vector<-rnorm(sim_n)
  S_T<-S*exp((r-0.5*sig^2)*maturity+sig*sqrt(maturity)*norm_vector)
  C_series<-exp(-r*maturity)*pmax(S_T-K,0)
  P_series<-exp(-r*maturity)*pmax(K-S_T,0)
  
  # Price
  Call<-mean(C_series)
  C_ste<-sd(C_series)/sqrt(sim_n)
  Put<-mean(P_series)
  P_ste<-sd(P_series)/sqrt(sim_n)
  
  return(list(Call,Put,C_ste,P_ste))
}

make_bs_mc_price_df<-function(S,K,r,maturity,sig,sim_n){
  price_list<-bs_mc_price(S,K,r,maturity,sig,sim_n)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  C_ste<-price_list[[3]]
  P_ste<-price_list[[4]]
  
  #Create dataframe
  mc_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(mc_cp_df)<-c('Call Price','Put Price')
  mc_ste_df<-data.frame('Call Ste'=C_ste,'Put Ste'=P_ste)
  colnames(mc_ste_df)<-c('Call Ste','Put Ste')
  
  return(list(mc_cp_df,mc_ste_df))
}


plot_bs_bt_euro_cov<-function(S,K,r,maturity,sig,tree_n){
  tree_n_vector<-seq(tree_n)
  bs_call_price<-bs_price(S,K,r,maturity,sig)[[1]]
  bs_put_price<-bs_price(S,K,r,maturity,sig)[[2]]
  bs_call_vector<-rep(bs_call_price,time=tree_n)
  bs_put_vector<-rep(bs_put_price,time=tree_n)
  bt_call_vector<-c()
  bt_put_vector<-c()
  for(n in tree_n_vector){
    bt_call_vector[n]<-bs_bt_euro_price(S,K,r,maturity,sig,n)[[1]]
    bt_put_vector[n]<-bs_bt_euro_price(S,K,r,maturity,sig,n)[[2]]
  }
  
  call_df1<-data.frame(Times=tree_n_vector,Call=bs_call_vector)
  call_df2<-data.frame(Times=tree_n_vector,Call=bt_call_vector)
  put_df1<-data.frame(Times=tree_n_vector,Put=bs_put_vector)
  put_df2<-data.frame(Times=tree_n_vector,Put=bt_put_vector)
  
  
  call_df<-call_df1%>% mutate(Type='BS') %>%
    bind_rows(call_df2%>%mutate(Type='BT'))
  
  put_df<-put_df1%>% mutate(Type='BS') %>%
    bind_rows(put_df2%>%mutate(Type='BT'))
  
  
  call_plot<-ggplot(data=call_df,aes(x=Times,y=Call,color=Type))+ 
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Scholes Model','Binomial Tree Model'))+
    labs(title='European Call Price : Black Scholes Model vs Binomial Tree Model',
         x='N Period',
         y='Call Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  
  put_plot<-ggplot(data=put_df,aes(x=Times,y=Put,color=Type))+
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Scholes Model','Binomial Tree Model'))+
    labs(title='European Put Price : Black Scholes Model vs Binomial Tree Model',
         x='N Period',
         y='Put Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  grid.arrange(call_plot,put_plot,nrow=1)
  
}

plot_bs_mc_euro_cov<-function(S,K,r,maturity,sig,sim_n){
  sim_n_vector<-seq(100,sim_n,100)
  bs_call_price<-bs_price(S,K,r,maturity,sig)[[1]]
  bs_put_price<-bs_price(S,K,r,maturity,sig)[[2]]
  bs_call_vector<-rep(bs_call_price,time=length(sim_n_vector))
  bs_put_vector<-rep(bs_put_price,time=length(sim_n_vector))
  
  mc_call_vector<-c()
  mc_put_vector<-c()
  
  
  for(i in seq(sim_n_vector)){
    mc_call_vector[i]<-bs_mc_price(S,K,r,maturity,sig,sim_n_vector[i])[[1]]
    mc_put_vector[i]<-bs_mc_price(S,K,r,maturity,sig,sim_n_vector[i])[[2]]
  }
  
  call_df1<-data.frame(Times=sim_n_vector,Call=bs_call_vector)
  call_df2<-data.frame(Times=sim_n_vector,Call=mc_call_vector)
  put_df1<-data.frame(Times=sim_n_vector,Put=bs_put_vector)
  put_df2<-data.frame(Times=sim_n_vector,Put=mc_put_vector)
  
  
  call_df<-call_df1%>% mutate(Type='BS') %>%
    bind_rows(call_df2%>%mutate(Type='MC'))
  
  put_df<-put_df1%>% mutate(Type='BS') %>%
    bind_rows(put_df2%>%mutate(Type='MC'))
  
  
  call_plot<-ggplot(data=call_df,aes(x=Times,y=Call,color=Type))+ 
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Scholes Model','Monte Carlo Simulation'))+
    labs(title='European Call Price : Black Scholes Model vs Monte Carlo Simulation',
         x='N Times',
         y='Call Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  
  put_plot<-ggplot(data=put_df,aes(x=Times,y=Put,color=Type))+
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Scholes Model','Monte Carlo Simulation'))+
    labs(title='European Put Price : Black Scholes Model vs Monte Carlo Simulation',
         x='N Times',
         y='Put Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  grid.arrange(call_plot,put_plot,nrow=1)
  
}


# tab2
bsm_price<-function(S,K,r,q,maturity,sig){
  if(maturity>0){
    d1<-(log(S/K)+(r-q+0.5*sig^2)*maturity)/(sig*sqrt(maturity))
    d2<-d1-sig*sqrt(maturity)
    N_d1<-pnorm(d1)
    n_d1<-dnorm(d1)
    N_d2<-pnorm(d2)
    n_d2<-dnorm(d2)
  }
  else{
    if(S>=K){
      N_d1<-1
      n_d1<-0
    }
    else{
      N_d1<-0
      n_d1<-0
    }
    N_d2<-N_d1
    n_d2<-n_d1
  }
  
  # Price
  Call<-S*exp(-q*maturity)*N_d1-K*exp(-r*maturity)*N_d2 
  Put<-K*exp(-r*maturity)*(1-N_d2)-S*exp(-q*maturity)*(1-N_d1)
  
  return(list(Call,Put,N_d1,N_d2,n_d1,n_d2))
}

make_bsm_price_df<-function(S,K,r,q,maturity,sig){
  price_list<-bsm_price(S,K,r,q,maturity,sig)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  N_d1<-price_list[[3]]
  N_d2<-price_list[[4]]
  n_d1<-price_list[[5]]
  n_d2<-price_list[[6]]
  
  N_md1<-1-N_d1
  N_md2<-1-N_d2
  
  # Greek
  delta_C<-exp(-q*maturity)*N_d1
  gamma_C<-exp(-q*maturity)*n_d1/(S*sig*sqrt(maturity))
  vega_C<-S*exp(-q*maturity)*n_d1*sqrt(T)
  rho_C<-K*maturity*exp(-r*maturity)*N_d2
  theta_C<--(S*exp(-q*maturity)*n_d1*sig/(2*sqrt(T))+K*r*exp(-r*maturity)*N_d2-q*S*exp(-q*maturity)*N_d1)
  
  delta_P<--exp(-q*maturity)*N_md1
  gamma_P<-gamma_C
  vega_P<-vega_C
  rho_P<--K*maturity*exp(-r*maturity)*N_md2
  theta_P<--(S*exp(-q*maturity)*n_d1*sig/(2*sqrt(T))-K*r*exp(-r*maturity)*N_md2+q*S*exp(-q*maturity)*N_md1)
  
  # Create dataframe
  bsm_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(bsm_cp_df)<-c('Call Price','Put Price')
  bsm_N_d_df<-data.frame('N(d1)'=N_d1,'N(d2)'=N_d2)
  colnames(bsm_N_d_df)<-c('N (d1)','N (d2)')
  bsm_greek_df<-data.frame(Delta=c(delta_C,delta_P),Gamma=c(gamma_C,gamma_P),Vega=c(vega_C,vega_P),
                          Rho=c(rho_C,rho_P),Theta=c(theta_C,theta_P))
  rownames(bsm_greek_df)<-c('Call','Put')
  
  return(list(bsm_cp_df,bsm_N_d_df,bsm_greek_df))
}


# Binomial Tree Model (European)
bsm_bt_euro_price<-function(S,K,r,q,maturity,sig,tree_n){
  t<-maturity/tree_n
  u<-exp(sig*sqrt(t))
  d<-1/u
  p<-(exp((r-q)*t)-d)/(u-d)
  q<-1-p
  FV_C<-0
  FV_P<-0
  
  for(i in seq(tree_n+1)){
    FV_C<-FV_C+choose(tree_n,i-1)*((p^(tree_n-i+1))*(q^(i-1)))*pmax(S*((u^(tree_n-i+1))*(d^(i-1)))-K,0)
    FV_P<-FV_P+choose(tree_n,i-1)*((p^(tree_n-i+1))*(q^(i-1)))*pmax(K-S*((u^(tree_n-i+1))*(d^(i-1))),0)
  }
  
  # Price
  Call<-exp(-r*maturity)*FV_C
  Put<-exp(-r*maturity)*FV_P
  
  return(list(Call,Put))
}

make_bsm_bt_euro_price_df<-function(S,K,r,q,maturity,sig,tree_n){
  price_list<-bsm_bt_euro_price(S,K,r,q,maturity,sig,tree_n)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  
  # Create dataframe
  bt_euro_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(bt_euro_cp_df)<-c('Call Price','Put Price')
  
  return(bt_euro_cp_df)
}


# Binomial Tree Model (American)
bsm_bt_amer_price<-function(S,K,r,q,maturity,sig,tree_n){
  t<-maturity/tree_n
  u<-exp(sig*sqrt(t))
  d<-1/u
  p<-(exp((r-q)*t)-d)/(u-d)
  q<-1-p
  discount<-exp(-r*t)
  
  lattice_C<-matrix(0L,nrow=tree_n+1,ncol=tree_n+1) 
  lattice_C[1,1]<-S
  
  # i:第幾層
  for(i in seq(tree_n)){
    # j:在i層中第幾個
    for(j in seq(i)){
      lattice_C[i+1,j+1]<-u*lattice_C[i,j]
    }
    lattice_C[i+1,1]<-d*lattice_C[i,1]
  }
  
  lattice_P<-lattice_C
  
  for(i in seq(tree_n+1,2)){
    for(j in seq(i,2)){
      if(i==tree_n+1){
        euro_C<-discount*(p*pmax(lattice_C[i,j]-K,0)+q*pmax(lattice_C[i,j-1]-K,0))
        euro_P<-discount*(p*pmax(K-lattice_P[i,j],0)+q*pmax(K-lattice_P[i,j-1],0))
      }
      else{
        euro_C<-discount*(p*lattice_C[i,j]+q*lattice_C[i,j-1])
        euro_P<-discount*(p*lattice_P[i,j]+q*lattice_P[i,j-1])
      }
      exer_C<-pmax(lattice_C[i-1,j-1]-K,0)
      exer_P<-pmax(K-lattice_P[i-1,j-1],0)
      lattice_C[i-1,j-1]<-pmax(euro_C,exer_C)
      lattice_P[i-1,j-1]<-pmax(euro_P,exer_P)
    }
  }
  
  # Price
  Call<-lattice_C[1,1]
  Put<-lattice_P[1,1]
  
  return(list(Call,Put,u,d,p))
  
}

make_bsm_bt_amer_price_df<-function(S,K,r,q,maturity,sig,tree_n){
  price_list<-bsm_bt_amer_price(S,K,r,q,maturity,sig,tree_n)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  u<-price_list[[3]]
  d<-price_list[[4]]
  p<-price_list[[5]]
  
  # Create dataframe
  bt_amer_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(bt_amer_cp_df)<-c('Call Price','Put Price')
  bt_arg_df<-data.frame('Up'=u,'Down'=d,'Up Probability'=p)
  colnames(bt_arg_df)<-c('Up','Down','Up Probability')
  
  return(list(bt_amer_cp_df,bt_arg_df))
}


# Monte Carlo Simulation
bsm_mc_price<-function(S,K,r,q,maturity,sig,sim_n){
  norm_vector<-rnorm(sim_n)
  S_T<-S*exp((r-q-0.5*sig^2)*maturity+sig*sqrt(maturity)*norm_vector)
  C_series<-exp(-r*maturity)*pmax(S_T-K,0)
  P_series<-exp(-r*maturity)*pmax(K-S_T,0)
  
  # Price
  Call<-mean(C_series)
  C_ste<-sd(C_series)/sqrt(sim_n)
  Put<-mean(P_series)
  P_ste<-sd(P_series)/sqrt(sim_n)
  
  return(list(Call,Put,C_ste,P_ste))
}

make_bsm_mc_price_df<-function(S,K,r,q,maturity,sig,sim_n){
  price_list<-bsm_mc_price(S,K,r,q,maturity,sig,sim_n)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  C_ste<-price_list[[3]]
  P_ste<-price_list[[4]]
  
  #Create dataframe
  mc_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(mc_cp_df)<-c('Call Price','Put Price')
  mc_ste_df<-data.frame('Call Ste'=C_ste,'Put Ste'=P_ste)
  colnames(mc_ste_df)<-c('Call Ste','Put Ste')
  
  return(list(mc_cp_df,mc_ste_df))
}


plot_bsm_bt_euro_cov<-function(S,K,r,q,maturity,sig,tree_n){
  tree_n_vector<-seq(tree_n)
  bs_call_price<-bsm_price(S,K,r,q,maturity,sig)[[1]]
  bs_put_price<-bsm_price(S,K,r,q,maturity,sig)[[2]]
  bs_call_vector<-rep(bs_call_price,time=tree_n)
  bs_put_vector<-rep(bs_put_price,time=tree_n)
  bt_call_vector<-c()
  bt_put_vector<-c()
  for(n in tree_n_vector){
    bt_call_vector[n]<-bsm_bt_euro_price(S,K,r,q,maturity,sig,n)[[1]]
    bt_put_vector[n]<-bsm_bt_euro_price(S,K,r,q,maturity,sig,n)[[2]]
  }
  
  call_df1<-data.frame(Times=tree_n_vector,Call=bs_call_vector)
  call_df2<-data.frame(Times=tree_n_vector,Call=bt_call_vector)
  put_df1<-data.frame(Times=tree_n_vector,Put=bs_put_vector)
  put_df2<-data.frame(Times=tree_n_vector,Put=bt_put_vector)
  
  
  call_df<-call_df1%>% mutate(Type='BS') %>%
    bind_rows(call_df2%>%mutate(Type='BT'))
  
  put_df<-put_df1%>% mutate(Type='BS') %>%
    bind_rows(put_df2%>%mutate(Type='BT'))
  
  
  call_plot<-ggplot(data=call_df,aes(x=Times,y=Call,color=Type))+ 
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Scholes Metron Model','Binomial Tree Model'))+
    labs(title='European Call Price : Black Scholes Metron Model vs Binomial Tree Model',
         x='N Period',
         y='Call Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  
  put_plot<-ggplot(data=put_df,aes(x=Times,y=Put,color=Type))+
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Scholes Metron Model','Binomial Tree Model'))+
    labs(title='European Put Price : Black Scholes Metron Model vs Binomial Tree Model',
         x='N Period',
         y='Put Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  grid.arrange(call_plot,put_plot,nrow=1)
  
}

plot_bsm_mc_euro_cov<-function(S,K,r,q,maturity,sig,sim_n){
  sim_n_vector<-seq(100,sim_n,100)
  bs_call_price<-bsm_price(S,K,r,q,maturity,sig)[[1]]
  bs_put_price<-bsm_price(S,K,r,q,maturity,sig)[[2]]
  bs_call_vector<-rep(bs_call_price,time=length(sim_n_vector))
  bs_put_vector<-rep(bs_put_price,time=length(sim_n_vector))
  
  mc_call_vector<-c()
  mc_put_vector<-c()
  
  
  for(i in seq(sim_n_vector)){
    mc_call_vector[i]<-bsm_mc_price(S,K,r,q,maturity,sig,sim_n_vector[i])[[1]]
    mc_put_vector[i]<-bsm_mc_price(S,K,r,q,maturity,sig,sim_n_vector[i])[[2]]
  }
  
  call_df1<-data.frame(Times=sim_n_vector,Call=bs_call_vector)
  call_df2<-data.frame(Times=sim_n_vector,Call=mc_call_vector)
  put_df1<-data.frame(Times=sim_n_vector,Put=bs_put_vector)
  put_df2<-data.frame(Times=sim_n_vector,Put=mc_put_vector)
  
  
  call_df<-call_df1%>% mutate(Type='BS') %>%
    bind_rows(call_df2%>%mutate(Type='MC'))
  
  put_df<-put_df1%>% mutate(Type='BS') %>%
    bind_rows(put_df2%>%mutate(Type='MC'))
  
  
  call_plot<-ggplot(data=call_df,aes(x=Times,y=Call,color=Type))+ 
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Scholes Metron Model','Monte Carlo Simulation'))+
    labs(title='European Call Price : Black Scholes Metron Model vs Monte Carlo Simulation',
         x='N Times',
         y='Call Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  
  put_plot<-ggplot(data=put_df,aes(x=Times,y=Put,color=Type))+
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Scholes Metron Model','Monte Carlo Simulation'))+
    labs(title='European Put Price : Black Scholes Metron Model vs Monte Carlo Simulation',
         x='N Times',
         y='Put Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  grid.arrange(call_plot,put_plot,nrow=1)
  
}


# tab3
# Black Model
bl_price<-function(Forward,K,r,maturity,sig){
  if(maturity>0){
    d1<-(log(Forward/K)+(0.5*sig^2)*maturity)/(sig*sqrt(maturity))
    d2<-d1-sig*sqrt(maturity)
    N_d1<-pnorm(d1)
    n_d1<-dnorm(d1)
    N_d2<-pnorm(d2)
    n_d2<-dnorm(d2)
  }
  else{
    if(Forward>=K){
      N_d1<-1
      n_d1<-0
    }
    else{
      N_d1<-0
      n_d1<-0
    }
    N_d2<-N_d1
    n_d2<-n_d1
  }
  
  # Price
  Call<-exp(-r*maturity)*(Forward*N_d1-K*N_d2) 
  Put<-exp(-r*maturity)*(K*(1-N_d2)-Forward*(1-N_d1))
  
  return(list(Call,Put,N_d1,N_d2,n_d1,n_d2))
}

make_bl_price_df<-function(Forward,K,r,maturity,sig){
  price_list<-bl_price(Forward,K,r,maturity,sig)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  N_d1<-price_list[[3]]
  N_d2<-price_list[[4]]
  n_d1<-price_list[[5]]
  n_d2<-price_list[[6]]
  
  N_md1<-1-N_d1
  N_md2<-1-N_d2
  
  # Greek
  delta_C<-exp(-r*maturity)*N_d1
  gamma_C<-exp(-r*maturity)*n_d1/(Forward*sig*sqrt(maturity))
  vega_C<-Forward*exp(-r*maturity)*n_d1*sqrt(T)
  rho_C<-K*maturity*exp(-r*maturity)*N_d2
  theta_C<--exp(-r*maturity)*(Forward*n_d1*sig/(2*sqrt(T))+K*r*N_d2-r*Forward*N_d1)
  
  delta_P<--exp(-r*maturity)*N_md1
  gamma_P<-gamma_C
  vega_P<-vega_C
  rho_P<--K*maturity*exp(-r*maturity)*N_md2
  theta_P<--exp(-r*maturity)*(Forward*n_d1*sig/(2*sqrt(T))-K*r*N_md2+r*Forward*N_md1)
  
  # Create dataframe
  bsm_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(bsm_cp_df)<-c('Call Price','Put Price')
  bsm_N_d_df<-data.frame('N(d1)'=N_d1,'N(d2)'=N_d2)
  colnames(bsm_N_d_df)<-c('N (d1)','N (d2)')
  bsm_greek_df<-data.frame(Delta=c(delta_C,delta_P),Gamma=c(gamma_C,gamma_P),Vega=c(vega_C,vega_P),
                           Rho=c(rho_C,rho_P),Theta=c(theta_C,theta_P))
  rownames(bsm_greek_df)<-c('Call','Put')
  
  return(list(bsm_cp_df,bsm_N_d_df,bsm_greek_df))
}


# Binomial Tree Model (European)
bl_bt_euro_price<-function(Forward,K,r,maturity,sig,tree_n){
  t<-maturity/tree_n
  u<-exp(sig*sqrt(t))
  d<-1/u
  p<-(1-d)/(u-d)
  q<-1-p
  FV_C<-0
  FV_P<-0
  
  for(i in seq(tree_n+1)){
    FV_C<-FV_C+choose(tree_n,i-1)*((p^(tree_n-i+1))*(q^(i-1)))*pmax(Forward*((u^(tree_n-i+1))*(d^(i-1)))-K,0)
    FV_P<-FV_P+choose(tree_n,i-1)*((p^(tree_n-i+1))*(q^(i-1)))*pmax(K-Forward*((u^(tree_n-i+1))*(d^(i-1))),0)
  }
  
  # Price
  Call<-exp(-r*maturity)*FV_C
  Put<-exp(-r*maturity)*FV_P
  
  return(list(Call,Put))
}

make_bl_bt_euro_price_df<-function(Forward,K,r,maturity,sig,tree_n){
  price_list<-bl_bt_euro_price(Forward,K,r,maturity,sig,tree_n)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  
  # Create dataframe
  bt_euro_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(bt_euro_cp_df)<-c('Call Price','Put Price')
  
  return(bt_euro_cp_df)
}


# Binomial Tree Model (American)
bl_bt_amer_price<-function(Forward,K,r,maturity,sig,tree_n){
  t<-maturity/tree_n
  u<-exp(sig*sqrt(t))
  d<-1/u
  p<-(1-d)/(u-d)
  q<-1-p
  discount<-exp(-r*t)
  
  lattice_C<-matrix(0L,nrow=tree_n+1,ncol=tree_n+1) 
  lattice_C[1,1]<-Forward
  
  # i:第幾層
  for(i in seq(tree_n)){
    # j:在i層中第幾個
    for(j in seq(i)){
      lattice_C[i+1,j+1]<-u*lattice_C[i,j]
    }
    lattice_C[i+1,1]<-d*lattice_C[i,1]
  }
  
  lattice_P<-lattice_C
  
  for(i in seq(tree_n+1,2)){
    for(j in seq(i,2)){
      if(i==tree_n+1){
        euro_C<-discount*(p*pmax(lattice_C[i,j]-K,0)+q*pmax(lattice_C[i,j-1]-K,0))
        euro_P<-discount*(p*pmax(K-lattice_P[i,j],0)+q*pmax(K-lattice_P[i,j-1],0))
      }
      else{
        euro_C<-discount*(p*lattice_C[i,j]+q*lattice_C[i,j-1])
        euro_P<-discount*(p*lattice_P[i,j]+q*lattice_P[i,j-1])
      }
      exer_C<-pmax(lattice_C[i-1,j-1]-K,0)
      exer_P<-pmax(K-lattice_P[i-1,j-1],0)
      lattice_C[i-1,j-1]<-pmax(euro_C,exer_C)
      lattice_P[i-1,j-1]<-pmax(euro_P,exer_P)
    }
  }
  
  # Price
  Call<-lattice_C[1,1]
  Put<-lattice_P[1,1]
  
  return(list(Call,Put,u,d,p))
  
}

make_bl_bt_amer_price_df<-function(Forward,K,r,maturity,sig,tree_n){
  price_list<-bl_bt_amer_price(Forward,K,r,maturity,sig,tree_n)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  u<-price_list[[3]]
  d<-price_list[[4]]
  p<-price_list[[5]]
  
  # Create dataframe
  bt_amer_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(bt_amer_cp_df)<-c('Call Price','Put Price')
  bt_arg_df<-data.frame('Up'=u,'Down'=d,'Up Probability'=p)
  colnames(bt_arg_df)<-c('Up','Down','Up Probability')
  
  return(list(bt_amer_cp_df,bt_arg_df))
}


# Monte Carlo Simulation
bl_mc_price<-function(Forward,K,r,maturity,sig,sim_n){
  norm_vector<-rnorm(sim_n)
  F_T<-Forward*exp((-0.5*sig^2)*maturity+sig*sqrt(maturity)*norm_vector)
  C_series<-exp(-r*maturity)*pmax(F_T-K,0)
  P_series<-exp(-r*maturity)*pmax(K-F_T,0)
  
  # Price
  Call<-mean(C_series)
  C_ste<-sd(C_series)/sqrt(sim_n)
  Put<-mean(P_series)
  P_ste<-sd(P_series)/sqrt(sim_n)
  
  return(list(Call,Put,C_ste,P_ste))
}

make_bl_mc_price_df<-function(Forward,K,r,maturity,sig,sim_n){
  price_list<-bl_mc_price(Forward,K,r,maturity,sig,sim_n)
  
  Call<-price_list[[1]]
  Put<-price_list[[2]]
  C_ste<-price_list[[3]]
  P_ste<-price_list[[4]]
  
  #Create dataframe
  mc_cp_df<-data.frame('Call Price'=Call,'Put Price'=Put)
  colnames(mc_cp_df)<-c('Call Price','Put Price')
  mc_ste_df<-data.frame('Call Ste'=C_ste,'Put Ste'=P_ste)
  colnames(mc_ste_df)<-c('Call Ste','Put Ste')
  
  return(list(mc_cp_df,mc_ste_df))
}


plot_bl_bt_euro_cov<-function(Forward,K,r,maturity,sig,tree_n){
  tree_n_vector<-seq(tree_n)
  bs_call_price<-bl_price(Forward,K,r,maturity,sig)[[1]]
  bs_put_price<-bl_price(Forward,K,r,maturity,sig)[[2]]
  bs_call_vector<-rep(bs_call_price,time=tree_n)
  bs_put_vector<-rep(bs_put_price,time=tree_n)
  bt_call_vector<-c()
  bt_put_vector<-c()
  for(n in tree_n_vector){
    bt_call_vector[n]<-bl_bt_euro_price(Forward,K,r,maturity,sig,n)[[1]]
    bt_put_vector[n]<-bl_bt_euro_price(Forward,K,r,maturity,sig,n)[[2]]
  }
  
  call_df1<-data.frame(Times=tree_n_vector,Call=bs_call_vector)
  call_df2<-data.frame(Times=tree_n_vector,Call=bt_call_vector)
  put_df1<-data.frame(Times=tree_n_vector,Put=bs_put_vector)
  put_df2<-data.frame(Times=tree_n_vector,Put=bt_put_vector)
  
  
  call_df<-call_df1%>% mutate(Type='BS') %>%
    bind_rows(call_df2%>%mutate(Type='BT'))
  
  put_df<-put_df1%>% mutate(Type='BS') %>%
    bind_rows(put_df2%>%mutate(Type='BT'))
  
  
  call_plot<-ggplot(data=call_df,aes(x=Times,y=Call,color=Type))+ 
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Model','Binomial Tree Model'))+
    labs(title='European Call Price : Black Model vs Binomial Tree Model',
         x='N Period',
         y='Call Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  
  put_plot<-ggplot(data=put_df,aes(x=Times,y=Put,color=Type))+
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Model','Binomial Tree Model'))+
    labs(title='European Put Price : Black Model vs Binomial Tree Model',
         x='N Period',
         y='Put Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  grid.arrange(call_plot,put_plot,nrow=1)
  
}

plot_bl_mc_euro_cov<-function(Forward,K,r,maturity,sig,sim_n){
  sim_n_vector<-seq(100,sim_n,100)
  bs_call_price<-bl_price(Forward,K,r,maturity,sig)[[1]]
  bs_put_price<-bl_price(Forward,K,r,maturity,sig)[[2]]
  bs_call_vector<-rep(bs_call_price,time=length(sim_n_vector))
  bs_put_vector<-rep(bs_put_price,time=length(sim_n_vector))
  
  mc_call_vector<-c()
  mc_put_vector<-c()
  
  
  for(i in seq(sim_n_vector)){
    mc_call_vector[i]<-bl_mc_price(Forward,K,r,maturity,sig,sim_n_vector[i])[[1]]
    mc_put_vector[i]<-bl_mc_price(Forward,K,r,maturity,sig,sim_n_vector[i])[[2]]
  }
  
  call_df1<-data.frame(Times=sim_n_vector,Call=bs_call_vector)
  call_df2<-data.frame(Times=sim_n_vector,Call=mc_call_vector)
  put_df1<-data.frame(Times=sim_n_vector,Put=bs_put_vector)
  put_df2<-data.frame(Times=sim_n_vector,Put=mc_put_vector)
  
  
  call_df<-call_df1%>% mutate(Type='BS') %>%
    bind_rows(call_df2%>%mutate(Type='MC'))
  
  put_df<-put_df1%>% mutate(Type='BS') %>%
    bind_rows(put_df2%>%mutate(Type='MC'))
  
  
  call_plot<-ggplot(data=call_df,aes(x=Times,y=Call,color=Type))+ 
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Model','Monte Carlo Simulation'))+
    labs(title='European Call Price : Black Model vs Monte Carlo Simulation',
         x='N Times',
         y='Call Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  
  put_plot<-ggplot(data=put_df,aes(x=Times,y=Put,color=Type))+
    geom_line()+
    scale_color_discrete(name='Type',labels=c('Black Model','Monte Carlo Simulation'))+
    labs(title='European Put Price : Black Model vs Monte Carlo Simulation',
         x='N Times',
         y='Put Price')+
    theme(plot.title=element_text(hjust=0.5,size=13))
  
  grid.arrange(call_plot,put_plot,nrow=1)
  
}








