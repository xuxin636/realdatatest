library(mvtnorm)
cond <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
E <- as.matrix(read.csv("/rigel/home/xx2319/realdatatest/0-1categoryofdata.csv"))
ww <- 40
w <- E[,2:(ww+1)]
J = ncol(w)
N = nrow(w)
K = 3

response <- w;
###my code###
Q <- matrix(1,J,K);Q[J,2:3] <- 0;Q[(J-1),3] <- 0;
##initial value###
A_initial <- matrix(0,J,K);A_initial[,1] <- runif(J,1,2);A_initial[,2] <- runif(J,1,2);A_initial[,3] <- runif(J,1,2);
A_initial <- A_initial*Q;
d_initial <- rnorm(J,0,1);D_initial <- cbind(d_initial,A_initial);
KK <- 20;theta_min <- -4;theta_max <- 4;mm1 <- seq(theta_min,theta_max,(theta_max-theta_min)/KK);mm <- mm1[-1]
THETA_tuta <- matrix(0,nrow=KK*KK*KK,ncol=3);THETA_tuta[,3] <- rep(mm,KK*KK);
THETA_tuta[,2] <-rep(c(rep(1,KK)%*%t(mm)),KK);THETA_tuta[,1] <-c(rep(1,KK*KK)%*%t(mm))#Õë¶ÔK <- 3µÄtheta·Ö¿é,»ñÈ¡thetaµÄ·Ö¿é
THETA_tuta <- cbind(rep(1,nrow(THETA_tuta)),THETA_tuta)
theta_square <- THETA_tuta[,2:4]*THETA_tuta[,2:4]
theta_tmp <- rowSums(theta_square)/2
xx <- seq(0.001,0.01,0.001);
xx1 <- matrix(0,nrow = length(xx)*length(xx)*length(xx),ncol=3);xx1[,3] <- rep(xx,length(xx)*length(xx));xx1[,2] <- rep(c(rep(1,length(xx))%*%t(xx)),length(xx))
xx1[,1] <- c(rep(1,length(xx)*length(xx))%*%t(xx))
lammda <- c(rep(0.01,20),rep(0.001,9),rep(0.001,10))*N;

soft <- function(a,b,K){
  for(k in 1:K){
    if(a[k]>0&a[k]>b[k]){a[k] <- a[k]-b[k]}
    else{a[k] <- 0}
    
  }
  return(a)
}
response <- t(response);
A_0 <- t(D_initial)
temp_0 <- like_temp_0 <- THETA_tuta%*%A_0
cc1 <- exp(like_temp_0%*%response-theta_tmp-rowSums(log(1+exp(like_temp_0))))
theta_post <- sweep(cc1, 2, colSums(cc1), "/") 
penalty<- -sum(sweep(A_0,2,lammda,"*"))
th0 <- rowSums(theta_post);
th1 <- THETA_tuta[,2]*th0;th2 <- THETA_tuta[,3]*th0;th3 <- THETA_tuta[,4]*th0;
uu <- theta_post%*%t(response);uu0 <-colSums(uu);uu1 <-colSums(THETA_tuta[,2]*uu);uu2 <-colSums(THETA_tuta[,3]*uu);uu3 <-colSums(THETA_tuta[,4]*uu);
frac_00 <- log(colSums(exp(temp_0%*%response-rowSums(log(1+exp(temp_0)))-theta_tmp)))
tol <- 1e-8
fan_0 <- sqrt(sum(A_0*A_0))
ss <- rep(0,502);ss[1] <- A_0[2,1]
timstart <- Sys.time()

xx <- 1/(1+exp(-temp_0))
A_grad <- uu0-colSums(th0* xx)
A_grad_2 <- -colSums(th0*xx*(1-xx))
d_tuta <- A_0[1,]-A_grad/A_grad_2
temp_1 <- THETA_tuta%*%rbind(d_tuta,A_0[-1,])
A_0[1,] <- d_tuta
for(k in 1:20){
  xx <- 1/(1+exp(-temp_1))
  A_grad <- uu1-colSums(th1* xx)
  A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,1])
  A_0[2,] <-A_0[2,]-A_grad/A_grad_2
  temp_1 <- THETA_tuta%*%A_0
}
xx <- 1/(1+exp(-temp_1))
A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,1])
A_0[2,]<-soft(A_0[2,],-lammda/A_grad_2,J)
temp_1 <- THETA_tuta%*%A_0

for(k in 1:20){
  xx <- 1/(1+exp(-temp_1))
  A_grad <- uu2-colSums(th2* xx)
  A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,2])
  A_0[3,1:(J-1)] <-A_0[3,1:(J-1)]-(A_grad/A_grad_2)[1:(J-1)];
  temp_1 <- THETA_tuta%*%A_0
}
xx <- 1/(1+exp(-temp_1))
A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,2])
temp_1 <- THETA_tuta%*%A_0
A_0[3,1:(J-1)] <- soft(A_0[3,1:(J-1)],(-lammda/A_grad_2)[1:(J-1)],J-1)

for(k in 1:20){
  xx <- 1/(1+exp(-temp_1))
  A_grad <- uu3-colSums(th3* xx)
  A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,3])
  A_0[4,1:(J-2)] <-A_0[4,1:(J-2)]-(A_grad/A_grad_2)[1:(J-2)];
  temp_1 <- THETA_tuta%*%A_0
}
xx <- 1/(1+exp(-temp_1))
A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,3])
A_0[4,1:(J-2)] <- soft(A_0[4,1:(J-2)],(-lammda/A_grad_2)[1:(J-2)],J-2)

ss[2] <- A_0[2,1]










####while####
for(m in 1:500){
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
  for(k in 1:20){
    xx <- 1/(1+exp(-temp_1))
    A_grad <- uu1-colSums(th1* xx)
    A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,1])
    A_0[2,] <-A_0[2,]-A_grad/A_grad_2
    temp_1 <- THETA_tuta%*%A_0
  }
  xx <- 1/(1+exp(-temp_1))
  A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,1])
  A_0[2,]<-soft(A_0[2,],-lammda/A_grad_2,J)
  temp_1 <- THETA_tuta%*%A_0
  
  for(k in 1:20){
    xx <- 1/(1+exp(-temp_1))
    A_grad <- uu2-colSums(th2* xx)
    A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,2])
    A_0[3,1:(J-1)] <-A_0[3,1:(J-1)]-(A_grad/A_grad_2)[1:(J-1)];
    temp_1 <- THETA_tuta%*%A_0
  }
  xx <- 1/(1+exp(-temp_1))
  A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,2])
  temp_1 <- THETA_tuta%*%A_0
  A_0[3,1:(J-1)] <- soft(A_0[3,1:(J-1)],(-lammda/A_grad_2)[1:(J-1)],J-1)
  
  for(k in 1:20){
    xx <- 1/(1+exp(-temp_1))
    A_grad <- uu3-colSums(th3* xx)
    A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,3])
    A_0[4,1:(J-2)] <-A_0[4,1:(J-2)]-(A_grad/A_grad_2)[1:(J-2)];
    temp_1 <- THETA_tuta%*%A_0
  }
  xx <- 1/(1+exp(-temp_1))
  A_grad_2 <- -colSums(th0*xx*(1-xx)*theta_square[,3])
  A_0[4,1:(J-2)] <- soft(A_0[4,1:(J-2)],(-lammda/A_grad_2)[-((J-1):J)],J-2)
  ss[m+2] <- A_0[2,1]
}

timend <- Sys.time()
tt <- timend-timstart
tt
bic <- -2*sum(log(colSums(exp(temp_0%*%response-rowSums(log(1+exp(temp_0)))-theta_tmp))))+log(N)*(J*K)
bic1 <- -2*sum(log(colSums(exp(temp_0%*%response-rowSums(log(1+exp(temp_0)))-theta_tmp))))+log(N)*(J*K)+2*sum(lammda*t(A_0))
RESULT <- rbind(c(bic,0,0,0),c(bic1,0,0,0),t(A_0))
write.csv(RESULT, file ='lammda01001001.csv')


