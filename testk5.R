library(mvtnorm)
cond <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
E <- as.matrix(read.csv("/rigel/home/xx2319/realdatatest/0-1categoryofdata.csv"))
ww <- 40
w <- E[,2:(ww+1)]

J = ncol(w)
N = nrow(w)
K = 5

response <- w;


###my code###
Q <- matrix(1,J,K);Q[36,2] <- 0;Q[c(3,36),3] <- 0;
##initial value###
A_initial <- matrix(0,J,K);A_initial[,1] <- runif(J,1,2);A_initial[,2] <- runif(J,1,2);A_initial[,3] <- runif(J,1,2);
A_initial <- A_initial;
d_initial <-sort(rnorm(J,0,1))[rank(colMeans(response))];

D_initial <- cbind(d_initial,A_initial);
KK <- 10;theta_min <- -4;theta_max <- 4;mm1 <- seq(theta_min,theta_max,(theta_max-theta_min)/KK);mm <- mm1[-1]


THETA_tuta <- matrix(0,nrow=KK^K,ncol=K);
for(k in 1:K){THETA_tuta[,K-k+1] <- rep(c(rep(1,KK^(k-1))%*%t(mm)),KK^(K-k))}
THETA_tuta <- cbind(rep(1,nrow(THETA_tuta)),THETA_tuta)
theta_square <- THETA_tuta[,2:(K+1)]*THETA_tuta[,2:(K+1)]
theta_tmp <- rowSums(theta_square)/2
xx <- seq(0.001,0.01,0.001);
xx1 <- matrix(0,nrow = length(xx)*length(xx)*length(xx),ncol=3);xx1[,3] <- rep(xx,length(xx)*length(xx));xx1[,2] <- rep(c(rep(1,length(xx))%*%t(xx)),length(xx))
xx1[,1] <- c(rep(1,length(xx)*length(xx))%*%t(xx))
lammda <- c(0.003,rep(0.003,18),0.003,0.009,rep(0.009,9),rep(0.0001,10))*N;
lammda <-rep(0,J)
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
cc1 <- exp(temp_0%*%response-theta_tmp-rowSums(log(1+exp(temp_0))))
theta_post <- sweep(cc1, 2, colSums(cc1), "/") 
th0 <- rowSums(theta_post);
th1 <- THETA_tuta[,2]*th0;th2 <- THETA_tuta[,3]*th0;th3 <- THETA_tuta[,4]*th0;
uu <- theta_post%*%t(response);uu0 <-colSums(uu);uu1 <-colSums(THETA_tuta[,2]*uu);uu2 <-colSums(THETA_tuta[,3]*uu);uu3 <-colSums(THETA_tuta[,4]*uu);
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
for(m in 1:300){
  temp_0 <- THETA_tuta%*%A_0
  cc1 <- exp(temp_0%*%response-theta_tmp-rowSums(log(1+exp(temp_0))))
  theta_post <- sweep(cc1, 2, colSums(cc1), "/") 
  th0 <- rowSums(theta_post);
  th1 <- THETA_tuta[,2]*th0;th2 <- THETA_tuta[,3]*th0;th3 <- THETA_tuta[,4]*th0;
  uu <- theta_post%*%t(response);uu0 <-colSums(uu);uu1 <-colSums(THETA_tuta[,2]*uu);uu2 <-colSums(THETA_tuta[,3]*uu);uu3 <-colSums(THETA_tuta[,4]*uu);
  
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
write.csv(RESULT, file =paste0('dim_j',cond,'.csv'))


