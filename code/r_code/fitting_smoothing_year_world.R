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
width = 238
tol = 0.008
disc_bound = 240
iter_tol = 1.0e-6

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


setwd("C:/Users/phans/Dropbox/temporal_awareness")


df <- read.csv(file="pooled_inter_year_world.csv", header=TRUE, sep=",") %>%
  as.data.frame() %>%
  filter(abs(L) <= width) 

#starting value
a <- 0.001
A <- -0.001
b <- 0.18
c <- 26
B <- -0.18

#Determine the numerator's polynomial order (best fit) 
numer_poly = 1
try(fit <- nls(beta ~ (b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L<=0 ) 
               , start=list(b= b, c = c, A = A, B = B), weights = w)
)
try(res1 <- sum(residuals(fit)^2))
if (exists("fit")==T){
  alter <- nls(beta ~ (a*L^2 + b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L<=0 ) 
               , start=list(a = a, b= b, c = c, A = A, B = B), weights = w)  
  res2 <- sum(residuals(alter)^2)
  dopoly2 <- try(res2<res1)
}
if (exists("fit")==F | dopoly2 ==T){
  numer_poly = 2
  fit <- nls(beta ~ (a*L^2 + b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L<=0 ) 
             , start=list(a = a, b= b, c = c, A = A, B = B), weights = w)  
}

pred <- predict(fit,df)
df_plot <- cbind(pred, df)

#compare with exponential
expo <- nls(beta ~ exp(Alpha + Beta*L), data =filter(df, L<=0), start=list(Alpha = 3, Beta= 0.08), weights = w)
pred_expo <- predict(expo,df)
df_plot <- cbind(pred_expo, df_plot)

#compare with hyperbolic
hype <- nls(beta ~ s/(1+k*L), data =filter(df, L<=0), start=list(s = 100, k= -0.15), weights = w)
pred_hype <- predict(hype,df)
df_plot <- cbind(pred_hype, df_plot)

#compare with polynomial
poly <- lm(beta ~poly(L, 8), data =filter(df, L<=0),
           weights = w)
pred_poly <- predict(poly,df)
df_plot <- cbind(pred_poly, df_plot)

#new starting value
coef <- coefficients(fit)[c("b","c","A","B")]
if (numer_poly == 2) {
  coef <- coefficients(fit)[c("a","b","c","A","B")]
}
x0<- coef 
DF <- filter(df, is.na(beta)==F & L <=0 )


if (numer_poly == 1) {
  optpar <- nloptr(x0=x0,
                   eval_f = sse,
                   eval_g_ineq=rootcons,
                   phase = "past",
                   L = DF$L,
                   y = DF$beta,
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
                             "print_level" = 1))
  optpar <- optpar$solution  
  rmse <- sum((DF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DF$L))^2) 
}

#iterate with new starting value
diff = 1
iteration = 0
while (diff >= tol) {
  iteration = iteration + 1
  x0 <- optpar
  
  if (numer_poly == 1) {
    optpar <- nloptr(x0=x0,
                     eval_f = sse,
                     eval_g_ineq=rootcons,
                     phase = "past",
                     L = DF$L,
                     y = DF$beta,
                     opts=list("algorithm"="NLOPT_LN_COBYLA",
                               maxeval=1e100,
                               "xtol_rel"=iter_tol,
                               "print_level" = 1))
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
                               "print_level" = 1))
    optpar <- optpar$solution  
    diff <- rmse - sum((DF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DF$L))^2)
    rmse <- sum((DF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DF$L))^2) 
    
  }  
  
}
print(paste0("RMSE is ",rmse," after ",iteration, " iterations"))

#project on discret value using best rational functional form
if (numer_poly == 1) {
  DF <- mutate(filter(df, L<=0), pred = (optpar[1]*L + optpar[2]) / (optpar[3] * L^2 + optpar[4] * L + 1))
  kernel_data <- mutate(kernel_data, pred = ifelse(L<=0,(optpar[1]*L + optpar[2]) / (optpar[3] * L^2 + optpar[4] * L + 1),0))
}
if (numer_poly == 2) {
  DF <- mutate(filter(df, L<=0), pred = (optpar[1]*L^2 + optpar[2]*L + optpar[3]) / (optpar[4] * L^2 + optpar[5] * L + 1))
  kernel_data <- mutate(kernel_data, pred = ifelse(L<=0,(optpar[1]*L^2 + optpar[2]*L + optpar[3]) / (optpar[4] * L^2 + optpar[5] * L + 1),0))
}
ggplot(df_plot) + geom_point(data = df_plot, aes(x = L, y = beta)) +
  geom_line(data = DF, aes(x = L, y = pred), col = "red") +
  geom_line(data = filter(df_plot, L<= 0 ), aes(x= L, y = pred_poly)) +
  geom_line(data = filter(df_plot, L<= 0 ), aes(x= L, y = pred_hype)) 


beta_p <- optpar
# futur ###########
a <- 0.00001
b <- -0.13
c <- 26.57
A <- -0.0006
B <- 0.11
x_basis <- c(a,b,c,A,B)
rm(fit)
dopoly2 = F
numer_poly = 1
try(fit <- nls(beta ~ (b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L>=0 )
               , start=list(b= b, c = c, A = A, B = B), weights = w)
)
try(res1 <- sum(residuals(fit)^2))
if (exists("fit")==T){
  alter <- nls(beta ~ (a*L^2 + b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L>=0 )
               , start=list(a = a, b= b, c = c, A = A, B = B), weights = w)
  res2 <- sum(residuals(alter)^2)
  dopoly2 <- try(res2<res1)
}

if (exists("fit")==F | dopoly2 ==T){
  numer_poly = 2
  fit <- nls(beta ~ (a*L^2 + b*L + c) / (A * L^2 + B * L + 1),data =filter(df, L>=0 )
             , start=list(a = a, b= b, c = c, A = A, B = B), weights = w)
}
pred <- predict(fit,df)
dff_plot <- cbind(pred, df)

#compare with exponential
expo <- nls(beta ~ exp(Alpha + Beta*L), data =filter(df, L>=0), start=list(Alpha = 3, Beta= -0.05), weights = w)
pred_expo <- predict(expo,df)
dff_plot <- cbind(pred_expo, dff_plot)

#compare with hyperbolic
hype <- nls(beta ~ s/(1+k*L), data =filter(df, L>=0), start=list(s = 100, k= 0.15), weights = w)
pred_hype <- predict(hype,df)
dff_plot <- cbind(pred_hype, dff_plot)

#compare with polynomial
poly <- lm(beta ~poly(L, 8), data =filter(df, L>=0),
           weights = w)
pred_poly <- predict(poly,df)
dff_plot <- cbind(pred_poly, dff_plot)


coef <- coefficients(fit)[c("b","c","A","B")]

if (numer_poly == 2) {
  coef <- coefficients(fit)[c("a","b","c","A","B")]
}

DFF <- filter(df, is.na(beta)==F & L >= 0 )


if (numer_poly == 1) {
  optpar <- nloptr(x0=x0,
                   eval_f = sse,
                   eval_g_ineq=rootcons,
                   phase = "futur",
                   L = DFF$L,
                   y = DFF$beta,
                   opts=list("algorithm"="NLOPT_LN_COBYLA",
                             maxeval=1e100,
                             "xtol_rel"=iter_tol,
                             "print_level" = 1))
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
                   opts=list("algorithm"="NLOPT_LN_COBYLA",
                             maxeval=1e100,
                             "xtol_rel"=iter_tol,
                             "print_level" = 1))
  optpar <- optpar$solution  
  rmse <- sum((DFF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DFF$L))^2) 
}


diff = 1
iteration = 0
while (diff >= tol) {
  iteration = iteration + 1
  x0 <- optpar
  
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
                               "print_level" = 1))
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
                     #ub = c(1, 10,100,-0.00000001,1),
                     opts=list("algorithm"="NLOPT_LN_COBYLA",
                               maxeval=1e100,
                               "xtol_rel"=iter_tol,
                               "print_level" = 1))
    optpar <- optpar$solution  
    diff <- rmse - sum((DFF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DFF$L))^2)
    rmse <- sum((DFF$beta-rational2(optpar[1],optpar[2],optpar[3],optpar[4],optpar[5],t=DFF$L))^2) 
    
  }  
  
}
print(paste0("RMSE is ",rmse," after ",iteration, " iterations"))

if (numer_poly == 1) {
  DFF <- mutate(filter(df, L>=0), pred = (optpar[1]*L + optpar[2]) / (optpar[3] * L^2 + optpar[4] * L + 1))
  kernel_data <- mutate(kernel_data, pred = ifelse(L>=0,(optpar[1]*L + optpar[2]) / (optpar[3] * L^2 + optpar[4] * L + 1), pred))
}
if (numer_poly == 2) {
  DFF <- mutate(filter(df, L>=0), pred = (optpar[1]*L^2 + optpar[2]*L + optpar[3]) / (optpar[4] * L^2 + optpar[5] * L + 1))
  kernel_data <- mutate(kernel_data, pred = ifelse(L>=0,(optpar[1]*L^2 + optpar[2]*L + optpar[3]) / (optpar[4] * L^2 + optpar[5] * L + 1),pred))
}
beta_f <- optpar
ggplot(df_plot) +
  geom_line(data = DF, aes(x = L, y = pred), color = "forestgreen", size = 1.5) +
  geom_line(data = DFF, aes(x = L, y = pred), color = "forestgreen", size = 1.5) +
  geom_line(data = filter(df_plot, L<=0), aes(x = L, y = pred), color = "forestgreen", size = 1, linetype = "dashed") +
  geom_line(data = filter(dff_plot, L>=0), aes(x = L, y = pred), color = "forestgreen", size = 1, linetype = "dashed") +
  geom_line(data = filter(df_plot, L<=0), aes(x = L, y = pred_expo), color = "red", size = 1) +
  geom_line(data = filter(dff_plot, L>=0), aes(x = L, y = pred_expo), color = "red", size = 1) +
  geom_vline(xintercept = 0, color = "maroon", linetype = "dashed", size = 1) +
  geom_hline(yintercept=0, colour="black", size=1) +
  geom_point(data = df_plot, aes(x = L, y = beta), color = "black", size = 3) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.minor.x = element_line(color = "grey80"),
    panel.grid.minor.y = element_line(color = "grey80"),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80")
  )


fitting_plot <- ggplot() +
  geom_line(data = filter(df_plot, L<=0), aes(x = L, y = pred_poly, color = "Polynomial"), size = 1, alpha = 0.6) +
  geom_line(data = filter(dff_plot, L>=0), aes(x = L, y = pred_poly, color = "Polynomial"), size = 1, alpha = 0.6) +   
  geom_line(data = filter(df_plot, L<=0), aes(x = L, y = pred_expo, color = "Exponential"), size = 1, alpha = 0.8) +
  geom_line(data = filter(dff_plot, L>=0), aes(x = L, y = pred_expo, color = "Exponential"), size = 1, alpha = 0.8) +
  geom_line(data = filter(df_plot, L<=0), aes(x = L, y = pred_hype, color = "Hyperbolic"), size = 1, alpha = 0.75) +
  geom_line(data = filter(dff_plot, L>=0), aes(x = L, y = pred_hype, color = "Hyperbolic"), size = 1, alpha = 0.75) +
  geom_line(data = DF, aes(x = L, y = pred, color = "Rational"), size = 1) +
  geom_line(data = DFF, aes(x = L, y = pred, color = "Rational"), size = 1) + 
  geom_point(data = df_plot, aes(x = L, y = beta), color = "black", size = 1.5,alpha = 0.85) +
  theme(
    legend.position = c(0, 1),legend.justification = c(-0.5, 1.5),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA)
  )  +
  labs(colour="",x="Days",y="Interest")+
  scale_color_manual(values = c("red", "deepskyblue4", "gold", "forestgreen"))

fitting_plot


ggsave(file = "fig3_kernel_yearly_world.pdf", plot = fitting_plot)

beta_matrix <- rbind(filter(df_plot, L<=0),filter(dff_plot, L>=0)) %>%
  select(-se, -min, -max, -w) %>%
  write.csv("fitted_kernel_yearly_world.csv")



kernel_data <- mutate(kernel_data, pred = pred / max(kernel_data$pred))
kern_plot <- ggplot(kernel_data) +
  geom_area(data = kernel_data, aes(x = L, y = pred), color ="black", fill = "black" ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA)
  )


kern_plot
ggsave(file = "kernel_inter_year_world.pdf", plot = kern_plot)

k_list = c()
lb_list = c()
ub_list = c()
interval_list = c()
intervals = c(-150,-56,-28,-7,0,7,28,56,150)
i = 1
while (i < 9) {

  expo <- nls(pred ~ exp(Alpha + Beta*L), data =filter(kernel_data, L<=intervals[i+1] & L>= intervals[i]), start=list(Alpha = 0, Beta= -0.08))
  # expo <- nls(pred ~ exp(Beta*L), data =filter(kernel_data, L<=inter[2] & L>= inter[1]), start=list( Beta= 0.08))
  # expo <- lm(pred ~ a*L+b, data = filter(kernel_data, L<=inter[2] & L>= inter[1]))

  pred_expo <- predict(expo,kernel_data)
  kernel_data_plot <- cbind(kernel_data, pred_expo)
  lb_list <- c(lb_list, confint2(expo)[2])
  ub_list <- c(ub_list, confint2(expo)[4])
  k_list <- c(k_list, coefficients(expo)[2])
  interval_list <- c(interval_list, paste0("[",intervals[i],",",intervals[i+1],"]"))
  i = i + 1
}
out_table <- interval_list %>%
  cbind(as.data.frame(k_list)) %>%
  cbind(lb_list) %>%
  cbind(ub_list)

colnames(out_table) <- c("time","k","lb","ub")
write.csv(out_table,"fig2_exponential_fit.csv", row.names=F)

