library(mvtnorm)
cond <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cond<-cond+2000
E <- as.matrix(read.csv("/rigel/home/xx2319/realdatatest/0-1categoryofdata.csv"))
ww <- 40
w <- E[,2:(ww+1)]
J = ncol(w)
N = nrow(w)
K = 3

response <- w;
###my code###
##initial value###
A_initial <- matrix(0,J,K);A_initial[,1] <- runif(J,1,2);A_initial[,2] <- runif(J,1,2);A_initial[,3] <- runif(J,1,2);
A_initial <- A_initial*Q;
d_initial <-sort(rnorm(J,0,1))[rank(colMeans(response))];
D_initial <- cbind(d_initial,A_initial);
KK <- 20;theta_min <- -4;theta_max <- 4;mm1 <- seq(theta_min,theta_max,(theta_max-theta_min)/KK);mm <- mm1[-1]
THETA_tuta <- matrix(0,nrow=KK*KK*KK,ncol=3);THETA_tuta[,3] <- rep(mm,KK*KK);
THETA_tuta[,2] <-rep(c(rep(1,KK)%*%t(mm)),KK);THETA_tuta[,1] <-c(rep(1,KK*KK)%*%t(mm))#Õë¶ÔK <- 3µÄtheta·Ö¿é,»ñÈ¡thetaµÄ·Ö¿é
THETA_tuta <- cbind(rep(1,nrow(THETA_tuta)),THETA_tuta)
theta_square <- THETA_tuta[,2:4]*THETA_tuta[,2:4]
theta_tmp <- rowSums(theta_square)/2
xx <- seq(0,0.03,0.0005);
xx1 <- matrix(0,nrow = length(xx)*length(xx),ncol=2);xx1[,2] <- rep(xx,length(xx));
xx1[,1] <- c(rep(1,length(xx))%*%t(xx))
lammda <- xx1[cond,]*N;

soft <- function(a,b,K){
  for(k in 1:K){
    if(a[k]>0&a[k]>b[k]){a[k] <- a[k]-b[k]}
    else{a[k] <- a[k]}
    
  }
  return(a)
}
response <- t(response);
A_0 <- t(D_initial)
temp_0  <- THETA_tuta%*%A_0
th0 <- 0;uu <- 0
for(i in 1:N){
  cc <- exp(rowSums(sweep(temp_0,2,response[,i]))-theta_tmp-rowSums(log(1+exp(temp_0))))
  th0 <-th0+cc/sum(cc)
  uu <- uu+(cc/sum(cc))%*%t(response[,i])
}

th1 <- THETA_tuta[,2]*th0;th2 <- THETA_tuta[,3]*th0;th3 <- THETA_tuta[,4]*th0;
uu0 <-colSums(uu);uu1 <-colSums(THETA_tuta[,2]*uu);uu2 <-colSums(THETA_tuta[,3]*uu);uu3 <-colSums(THETA_tuta[,4]*uu);
ss <- rep(0,502);ss[1] <- A_0[2,1]
timstart <- Sys.time()

xx <- 1/(1+exp(-temp_0))
A_grad <- uu0-colSums(th0* xx)
A_grad_2 <- -colSums(th0*xx*(1-xx))
d_tuta <- A_0[1,]-A_grad/A_grad_2
temp_1 <- THETA_tuta%*%rbind(d_tuta,A_0[-1,])
A_0[1,] <- d_tuta

for(k in 2:(K+1)){
  
  for(s in 1:20){
    xx <- 1/(1+exp(-temp_1))
    A_grad <- colSums(THETA_tuta[,k]*uu)-colSums(THETA_tuta[,k]*th0* xx)
    A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,(k-1)])
    A_0[k,] <-A_0[k,]-A_grad/A_grad_2
    temp_1 <- THETA_tuta%*%A_0
  }
  xx <- 1/(1+exp(-temp_1))
  A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,(k-1)])
  A_0[2,]<-soft(A_0[2,],-lammda/A_grad_2,J)
  temp_1 <- THETA_tuta%*%A_0
}


ss[2] <- A_0[2,1]

####while####
for(m in 1:500){
  temp_0  <- THETA_tuta%*%A_0
  th0 <- 0;uu <- 0
  for(i in 1:N){
    cc <- exp(rowSums(sweep(temp_0,2,response[,i]))-theta_tmp-rowSums(log(1+exp(temp_0))))
    th0 <-th0+cc/sum(cc)
    uu <- uu+(cc/sum(cc))%*%t(response[,i])
  }
  
  th1 <- THETA_tuta[,2]*th0;th2 <- THETA_tuta[,3]*th0;th3 <- THETA_tuta[,4]*th0;
  uu0 <-colSums(uu);uu1 <-colSums(THETA_tuta[,2]*uu);uu2 <-colSums(THETA_tuta[,3]*uu);uu3 <-colSums(THETA_tuta[,4]*uu);  
  xx <- 1/(1+exp(-temp_0))
  A_grad <- uu0-colSums(th0* xx)
  A_grad_2 <- -colSums(th0*xx*(1-xx))
  d_tuta <- A_0[1,]-A_grad/A_grad_2
  temp_1 <- THETA_tuta%*%rbind(d_tuta,A_0[-1,])
  A_0[1,] <- d_tuta
  
  for(k in 2:(K+1)){
    
    for(s in 1:20){
      xx <- 1/(1+exp(-temp_1))
      A_grad <- colSums(THETA_tuta[,k]*uu)-colSums(THETA_tuta[,k]*th0* xx)
      A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,(k-1)])
      A_0[k,] <-A_0[k,]-A_grad/A_grad_2
      temp_1 <- THETA_tuta%*%A_0
    }
    xx <- 1/(1+exp(-temp_1))
    A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,(k-1)])
    A_0[2,]<-soft(A_0[2,],-lammda/A_grad_2,J)
    temp_1 <- THETA_tuta%*%A_0
  }
  ss[m+2] <- A_0[2,1]
}

timend <- Sys.time()
tt <- timend-timstart
tt
bic <- -2*sum(log(colSums(exp(temp_0%*%response-rowSums(log(1+exp(temp_0)))-theta_tmp))))+log(N)*(J*K)
bic1 <- -2*sum(log(colSums(exp(temp_0%*%response-rowSums(log(1+exp(temp_0)))-theta_tmp))))+log(N)*(J*K)+2*sum(lammda*t(A_0))
RESULT <- rbind(c(bic,0,0,0),c(bic1,0,0,0),t(A_0))
write.csv(RESULT, file =paste0('dim_n',cond,'.csv'))
