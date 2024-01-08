library(coin)
library(snow)
library(tidyverse)
library(MKinfer)
library(future)
library(future.apply)



results <- read.csv("fullsamples_test.csv") %>% 
  mutate(sp = as.factor(sp))

permutation_test100 <- function(results, n ) {
  p_frame <- data.frame(matrix(nrow = 1, ncol = 2))
  names(p_frame) <- c("p", "sp")
  p_data <- results2 %>% 
    filter(sp == n | sp %in% seq(n, n + 100,1)) %>% 
    mutate(sp = ifelse(!sp == n, n+100, sp)) %>% 
    dplyr::select(cor.sp, sp) 
  ttest <- perm.t.test(cor.sp ~ sp, data = p_data,  R = 999, alternative = "less")
  p_frame$perm_p <- ttest$perm.p.value
  p_frame$p <-ttest$p.value
  p_frame$sp <- paste(unique(p_data$sp)[1],"-",unique(p_data$sp[2]))
  return(p_frame)
}

permutation_test <- function(results, n ) {
  p_frame <- data.frame(matrix(nrow = 1, ncol = 2))
  names(p_frame) <- c("p", "sp")
  p_data <- results2 %>% 
    filter(sp == n | sp == n + 100) %>% 
    #mutate(sp = ifelse(!sp == n, n+100, sp)) %>% 
    dplyr::select(cor.sp, sp) 
  ttest <- perm.t.test(cor.sp ~ sp, data = p_data,  R = 999, alternative = "less")
  p_frame$perm_p <- ttest$perm.p.value
  p_frame$p <-ttest$p.value
  p_frame$sp <- paste(unique(p_data$sp)[1],"-",unique(p_data$sp[2]))
  return(p_frame$perm_p)
}

#as.data.table(results2)[, perm.t.test(cor.sp ~ sp, data = p_data,  R = 999, alternative = "less"), by = c("id", "sp")]


resultsxx <- results %>% 
  filter(id %in% c("3")) %>% 
  filter(sp == 100)




permutation_test_cont <- function(results, n ) {

  p_data <- results %>% 
    group_by(id) %>% 
    filter(sp == n | sp == n + 100) %>% # sp == ceiling(n + n*0.25)
    dplyr::select(cor.sp, sp, id) 
  
  unique_sp <- unique(p_data$sp)
  
  test_data <- p_data %>%
    mutate(steps = paste(unique_sp[1], "-", unique_sp[2]),
           sp = as.factor(sp)) %>%
    group_by(id) %>%
    reframe(perm_p = coin::pvalue(oneway_test(cor.sp ~ sp,  
                                        distribution = approximate(nresample = 9999),
                                                                   #parallel = "snow",
                                                                   #cl = cl), 
                                        alternative = "less"))[1], 
            norm_p = t.test(cor.sp ~ sp,  alternative = "less")$p.value,
            n = n,
            steps = unique(steps), 
            id = unique(id)) 
  
  return(test_data)
}

t_test_cont <- function(results, n ) {
  
  p_data <- results %>% 
    group_by(id) %>% 
    filter(sp == n | sp == n + 100) %>% 
    #mutate(sp = ifelse(!sp == n, n+100, sp)) %>% 
    # sp == ceiling(n + n*0.25)
    dplyr::select(cor.sp, sp, id) 
  
  unique_sp <- unique(p_data$sp)
  
  test_data <- p_data %>%
    mutate(steps = paste(unique_sp[1], "-", unique_sp[2]),
           sp = as.factor(sp)) %>%
    group_by(id) %>%
    reframe(norm_p = t.test(cor.sp ~ sp, alternative = "less")$p.value,
            n = n,
            steps = unique(steps), 
            id = unique(id)) 
  return(test_data)
}

plan(multisession, workers = 6)

start <- Sys.time()
a <- future_lapply(seq(100,2500, 1), FUN = permutation_test_cont, results = results, future.seed = T)
#b <- lapply(1:1500, FUN = permutation_test100, results = results)
Sys.time() - start

xx <- do.call(rbind,a)


plot(xx$norm_p ~ xx$n)
plot(xx$norm_p ~ xx$n)
lines(x=c(0,1900), y= c(0.05, 0.05), lwd =3, col = "red")
boxplot(xx$norm_p ~ xx$id)
hist(xx$norm_p)
lmm <- lm(norm_p ~ n, data = xx)
summary(lmm)
plot(lmm)
test <- xx %>% 
  group_by(id) %>% 
  filter(perm_p >= 0.05) %>% 
  slice_min(n = 1, order_by = n) 

aa <- unique(test$n)
plot(aa)
boxplot(aa)

test2 <- results %>% 
  filter(sp %in% c(test$n, (test$n * 0.25))) %>% 
  group_by(sp, id) %>% 
  reframe(mean = mean(cor.sp, na.rm = T))
