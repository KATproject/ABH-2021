library(pacman)
pacman::p_load(foreign
               ,reshape2
               ,tidyr
               ,dplyr
               ,ggplot2
               ,stringr
               ,lubridate
               ,forecast
               ,lfe
               ,ggpubr
               ,officer
               ,directlabels
               ,nloptr
               ,nlstools

)
rm(list=ls(all=TRUE))

#options(scipen = 999)
print_detail = 0
width = 238
tol = 0.008
dopoly2 = F
disc_bound = 240
iter_tol = 1.0e-6

bs_min = 1
bs_max = 2000

kernel_data <- data.frame(L = -width:width) 

### Optimization program  ### ----
rational=function(b,c,A,B,t){
  prd=(b*t + c) / (A*t^2 + B*t + 1)
  return(prd)
}
rational2=function(a,b,c,A,B,t){
  prd=(a*t^2 + b*t + c) / (A*t*t + B*t + 1)
  return(prd)
}
sse <- function(gradient,L,y,phase){
  sum((y - rational(gradient[1],gradient[2],gradient[3],gradient[4],t=L))^2)
}
sse2 <- function(gradient,L,y,phase){
  sum((y - rational2(gradient[1],gradient[2],gradient[3],gradient[4],gradient[5],t=L))^2)
}

rootcons <- function(gradient,L,y,phase) {
  delta <- gradient[4]^2 - (4*gradient[3])
  if (delta > 0 & gradient[3] != 0) {
    x1 <- (-gradient[4] - sqrt(delta)) / (2*gradient[3])
    x2 <- (-gradient[4] + sqrt(delta)) / (2*gradient[3])
    
    if ((x1>0 & x2>0 & phase == "past") | (x1<0 & x2<0 & phase == "futur")) {
      roots <- 0
    }
    
    if ((x1<0 & x2<0 & phase == "past") | (x1>0 & x2>0 & phase == "futur")) {
      roots <- numeric(2)
      if (phase == "past") {
        roots[1] <- disc_bound + ((-gradient[4] + sqrt(gradient[4]^2 - (4*gradient[3]))) / (2*gradient[3])) 
        roots[2] <- disc_bound + ((-gradient[4] - sqrt(gradient[4]^2 - (4*gradient[3]))) / (2*gradient[3]))           
      } else {
        roots[1] <- disc_bound - ((-gradient[4] + sqrt(gradient[4]^2 - (4*gradient[3]))) / (2*gradient[3])) 
        roots[2] <- disc_bound - ((-gradient[4] - sqrt(gradient[4]^2 - (4*gradient[3]))) / (2*gradient[3]))           
      }
    }
    
    if ((x1<0 & x2>0) | (x1>0 & x2<0)) {
      roots <- numeric(2)
      X <- c(x1,x2)
      neg <- which.min(X)
      pos <- which.max(X)
      if (phase == "past") {
        roots[1] <- disc_bound + X[neg]
        roots[2] <- -X[pos]        
      } else {
        roots[1] <- X[neg]
        roots[2] <- disc_bound - X[pos]          
      }
    }      
  } else{
    roots <- 0
  }
  return(roots)
}

rootcons2 <- function(gradient,L,y,phase) {
  delta <- gradient[5]^2 - (4*gradient[4])
  if (delta > 0 & gradient[4] != 0) {
    x1 <- (-gradient[5] - sqrt(delta)) / (2*gradient[4])
    x2 <- (-gradient[5] + sqrt(delta)) / (2*gradient[4])
    
    if ((x1>0 & x2>0 & phase == "past") | (x1<0 & x2<0 & phase == "futur")) {
      roots <- 0
    }
    
    if ((x1<0 & x2<0 & phase == "past") | (x1>0 & x2>0 & phase == "futur")) {
      roots <- numeric(2)
      if (phase == "past") {
        roots[1] <- disc_bound + ((-gradient[5] + sqrt(gradient[5]^2 - (4*gradient[4]))) / (2*gradient[4])) 
        roots[2] <- disc_bound + ((-gradient[5] - sqrt(gradient[5]^2 - (4*gradient[4]))) / (2*gradient[4]))           
      } else {
        roots[1] <- disc_bound - ((-gradient[5] + sqrt(gradient[5]^2 - (4*gradient[4]))) / (2*gradient[4])) 
        roots[2] <- disc_bound - ((-gradient[5] - sqrt(gradient[5]^2 - (4*gradient[4]))) / (2*gradient[4]))           
      }
    }
    
    if ((x1<0 & x2>0) | (x1>0 & x2<0)) {
      roots <- numeric(2)
      X <- c(x1,x2)
      neg <- which.min(X)
      pos <- which.max(X)
      if (phase == "past") {
        roots[1] <- disc_bound + X[neg]
        roots[2] <- -X[pos]        
      } else {
        roots[1] <- X[neg]
        roots[2] <- disc_bound - X[pos]          
      }
    }      
  } else{
    roots <- 0
  }
  return(roots)
}


############################


try(setwd("C:/Users/phans/Dropbox/temporal_awareness"))
try(setwd("/home/sphan/Research_projects/temporal_awareness"))



df <- read.csv(file="saves/BS_kat_1234.csv", header=TRUE, sep=",") %>%
  as.data.frame() %>% 
  mutate(L=L*7) %>% 
  mutate(w = 1/se) %>%
  filter(abs(L) <= width) %>%
  filter(ite >= bs_min & ite <= bs_max)
# opt_p = c(-0.003112885, -0.611805995, 20.784559678, -0.002073140, -0.744111490)
# opt_f = c(-0.0008652961,  0.1079526320, 20.7371270573, -0.0020379253, 0.5109108137)

# a <- opt_p[1]
# b <- opt_p[2]
# c <- opt_p[3]
# A <- opt_p[4]
# B <- opt_p[5]
a <- 0.001
A <- -0.001
b <- 0.18
c <- 26
B <- -0.18

speci = bs_min + 1000

numer_poly = 1
try(fit <- nls(beta ~ (b*L + c) / (A * L^2 + B * L + 1), data =filter(df, L<=0 & ite==speci)
               , start=list(b= b, c = c, A = A, B = B), weights = w)
)
try(res1 <- sum(residuals(fit)^2))
if (exists("fit")==T){
  alter <- nls(beta ~ (a*L^2 + b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L<=0 & ite==speci )
               , start=list(a = a, b= b, c = c, A = A, B = B), weights = w)
  res2 <- sum(residuals(alter)^2)
  dopoly2 <- try(res2<res1)
}
if (exists("fit")==F | dopoly2 ==T){
  numer_poly = 2
  fit <- nls(beta ~ (a*L^2 + b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L<=0 & ite==speci)
             , start=list(a = a, b= b, c = c, A = A, B = B), weights = w)
}

pred <- predict(fit,df)
df_plot <- cbind(pred, df)

coef <- coefficients(fit)[c("b","c","A","B")]
if (numer_poly == 2) {
  coef <- coefficients(fit)[c("a","b","c","A","B")]
}
x0<- coef
DF <- filter(df, is.na(beta)==F & L <=0 & ite==speci)




if (numer_poly == 1) {
  optpar <- nloptr(x0=x0,
                   eval_f = sse,
                   eval_g_ineq=rootcons,
                   phase = "past",
                   L = DF$L,
                   y = DF$beta,
                   # ub = c(10,100,-0.00000001,1),
                   opts=list("algorithm"="NLOPT_LN_COBYLA",
                             maxeval=1e100,
                             "xtol_rel"=iter_tol,
                             "print_level" = 1))
  optpar <- optpar$solution  
  rmse <- sum((DF$beta-rational(optpar[1],optpar[2],optpar[3],optpar[4],t=DF$L))^2) 
}

if (numer_poly == 2) {
  optpar <- nloptr(x0=x0,
                   eval_f = sse2,
                   eval_g_ineq=rootcons2,
                   phase = "past",
                   L = DF$L,
                   y = DF$beta,
                   ub = c(1, 10,100,-0.00000001,-0.0000001),
                   opts=list("algorithm"="NLOPT_LN_COBYLA",
                             maxeval=1e100,
                             "xtol_rel"=iter_tol,
                             "print_level" = print_detail))
  optpar <- optpar$solution  
  rmse <- sum((DF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DF$L))^2) 
}
optpar_base <- optpar
for (val in bs_min:bs_max) {
  print(paste0("iteration ",val))
  DF <- filter(df, is.na(beta)==F & L <=0 & ite == val )
  diff = 1
  iteration = 0
  while (diff >= tol) {
    iteration = iteration + 1
    if (iteration == 1) {
      x0 <- optpar_base
    }
    else {
      x0 <- optpar
    }
    if (numer_poly == 1) {
      optpar <- nloptr(x0=x0,
                       eval_f = sse,
                       eval_g_ineq=rootcons,
                       phase = "past",
                       L = DF$L,
                       y = DF$beta,
                       # ub = c(10,100,-0.00000001,-0.0000001),
                       opts=list("algorithm"="NLOPT_LN_COBYLA",
                                 maxeval=1e100,
                                 "xtol_rel"=iter_tol,
                                 "print_level" = print_detail))
      optpar <- optpar$solution  
      diff <- rmse - sum((DF$beta-rational(optpar[1],optpar[2],optpar[3],optpar[4],t=DF$L))^2)
      rmse <- sum((DF$beta-rational(optpar[1],optpar[2],optpar[3],optpar[4],t=DF$L))^2) 
    }
    
    if (numer_poly == 2) {
      optpar <- nloptr(x0=x0,
                       eval_f = sse2,
                       eval_g_ineq=rootcons2,
                       phase = "past",
                       L = DF$L,
                       y = DF$beta,
                       ub = c(1, 10,100,-0.00000001,-0.0000001),
                       opts=list("algorithm"="NLOPT_LN_COBYLA",
                                 maxeval=1e100,
                                 "xtol_rel"=1.0e-6,
                                 "print_level" = print_detail))
      optpar <- optpar$solution  
      diff <- rmse - sum((DF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DF$L))^2)
      rmse <- sum((DF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DF$L))^2) 
      
    }  
  } #END minimization loop
  assign(paste0("optpar",val),optpar)
} #END fixed effect loop

for (val in bs_min:bs_max) {
  print(paste0("pred ",val))
  optpar <- get(paste0("optpar",val))
  varN <- paste0("pred", val)
  if (numer_poly == 1) {
    DFtmp <- mutate(kernel_data, temppred = (optpar[1]*L + optpar[2]) / (optpar[3] * L^2 + optpar[4] * L + 1))
    kernel_data <- mutate(kernel_data, temppred = ifelse(L<=0,(optpar[1]*L + optpar[2]) / (optpar[3] * L^2 + optpar[4] * L + 1),0))
  }
  if (numer_poly == 2) {
    DFtmp <- mutate(kernel_data, temppred = (optpar[1]*L^2 + optpar[2]*L + optpar[3]) / (optpar[4] * L^2 + optpar[5] * L + 1))
    kernel_data <- mutate(kernel_data, temppred = ifelse(L<=0,(optpar[1]*L^2 + optpar[2]*L + optpar[3]) / (optpar[4] * L^2 + optpar[5] * L + 1),0))
  }
  kernel_data[[varN]] <- DFtmp$temppred
  rm(list = paste0("optpar",val))
}
rm(DFtmp)



# futur ###########
# a <- opt_f[1]
# b <- opt_f[2]
# c <- opt_f[3]
# A <- opt_f[4]
# B <- opt_f[5]

a <- 0.00001
b <- -0.13
c <- 26.57
A <- -0.0006
B <- 0.11
x_basis <- c(a,b,c,A,B)
rm(fit)
dopoly2 = F
numer_poly = 1
try(fit <- nls(beta ~ (b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L>=0 & ite==speci)
               , start=list(b= b, c = c, A = A, B = B), weights = w)
)
try(res1 <- sum(residuals(fit)^2))
if (exists("fit")==T){
  alter <- nls(beta ~ (a*L^2 + b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L>=0 & ite==speci)
               , start=list(a = a, b= b, c = c, A = A, B = B), weights = w)
  res2 <- sum(residuals(alter)^2)
  dopoly2 <- try(res2<res1)
}

if (exists("fit")==F | dopoly2 ==T){
  numer_poly = 2
  fit <- nls(beta ~ (a*L^2 + b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L>=0 & ite==speci)
             , start=list(a = a, b= b, c = c, A = A, B = B), weights = w)
}

coef <- coefficients(fit)[c("b","c","A","B")]

if (numer_poly == 2) {
  coef <- coefficients(fit)[c("a","b","c","A","B")]
}


x0<-coef
DFF <- filter(df, is.na(beta)==F & L >= 0 & ite==speci)


if (numer_poly == 1) {
  optpar <- nloptr(x0=x0,
                   eval_f = sse,
                   eval_g_ineq=rootcons,
                   phase = "futur",
                   L = DFF$L,
                   y = DFF$beta,
                   #ub = c(10,100,-0.00000001,1),
                   opts=list("algorithm"="NLOPT_LN_COBYLA",
                             maxeval=1e100,
                             "xtol_rel"=iter_tol,
                             "print_level" = print_detail))
  optpar <- optpar$solution  
  rmse <- sum((DFF$beta-rational(optpar[1],optpar[2],optpar[3],optpar[4],t=DFF$L))^2) 
}

if (numer_poly == 2) {
  optpar <- nloptr(x0=x0,
                   eval_f = sse2,
                   eval_g_ineq=rootcons2,
                   phase = "futur",
                   L = DFF$L,
                   y = DFF$beta,
                   #ub = c(1, 10,100,-0.00000001,1),
                   opts=list("algorithm"="NLOPT_LN_COBYLA",
                             maxeval=1e100,
                             "xtol_rel"=iter_tol,
                             "print_level" = print_detail))
  optpar <- optpar$solution  
  rmse <- sum((DFF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DFF$L))^2) 
}
optpar_base <- optpar
rm(optar)
for (val in bs_min:bs_max) {
  DFF <- filter(df, is.na(beta)==F & L >=0 & ite == val )
  diff = 1
  iteration = 0
  print(paste0("iteration (futur) ",val))
  while (diff >= tol) {
    iteration = iteration + 1
    if (iteration == 1) {
      x0 <- optpar_base
    }
    else {
      x0 <- optpar
    }
    if (numer_poly == 1) {
      optpar <- nloptr(x0=x0,
                       eval_f = sse,
                       eval_g_ineq=rootcons,
                       phase = "futur",
                       L = DFF$L,
                       y = DFF$beta,
                       #ub = c(10,100,-0.00000001,1),
                       opts=list("algorithm"="NLOPT_LN_COBYLA",
                                 maxeval=1e100,
                                 "xtol_rel"=iter_tol,
                                 "print_level" = print_detail))
      optpar <- optpar$solution  
      diff <- rmse - sum((DFF$beta-rational(optpar[1],optpar[2],optpar[3],optpar[4],t=DFF$L))^2)
      rmse <- sum((DFF$beta-rational(optpar[1],optpar[2],optpar[3],optpar[4],t=DFF$L))^2) 
    }
    
    if (numer_poly == 2) {
      optpar <- nloptr(x0=x0,
                       eval_f = sse2,
                       eval_g_ineq=rootcons2,
                       phase = "futur",
                       L = DFF$L,
                       y = DFF$beta,
                       ub = c(1, 10,100,-0.00000001,1),
                       opts=list("algorithm"="NLOPT_LN_COBYLA",
                                 maxeval=1e100,
                                 "xtol_rel"=iter_tol,
                                 "print_level" = print_detail))
      optpar <- optpar$solution  
      diff <- rmse - sum((DFF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DFF$L))^2)
      rmse <- sum((DFF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DFF$L))^2) 
      
    }  
  }
  assign(paste0("optpar",val),optpar)
}




for (val in bs_min:bs_max) {
  optpar <- get(paste0("optpar",val))
  varN <- paste0("pred", val)
  print(paste0("pred (futur) ",val))
  
  if (numer_poly == 1) {
    #DFtmp <- mutate(kernel_data, temppred = (optpar[1]*L + optpar[2]) / (optpar[3] * L^2 + optpar[4] * L + 1))
    kernel_data <- mutate(kernel_data, !!varN := ifelse(L>=0,(optpar[1]*L + optpar[2]) / (optpar[3] * L^2 + optpar[4] * L + 1),get(varN)))
  }
  if (numer_poly == 2) {
    #DFtmp <- mutate(kernel_data, temppred = (optpar[1]*L^2 + optpar[2]*L + optpar[3]) / (optpar[4] * L^2 + optpar[5] * L + 1))
    kernel_data <- mutate(kernel_data, !!varN := ifelse(L>=0,(optpar[1]*L^2 + optpar[2]*L + optpar[3]) / (optpar[4] * L^2 + optpar[5] * L + 1),get(varN)))
  }
  rm(list = paste0("optpar",val))
  #kernel_data[[varN]] <- DFtmp$temppred
}

beta_matrix <- rbind(filter(kernel_data, L<=0),filter(kernel_data, L>=0)) %>%
  select(-temppred) %>%
  write.csv(paste0("saves/fitted_kernel_daily_BS_seeds",bs_min,"_",bs_max,".csv"))
