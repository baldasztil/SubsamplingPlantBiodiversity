library(coin)
library(snow)
library(tcontinentyverse)
library(MKinfer)
library(future)
library(future.apply)
library(arrow)
library(data.table)



path <- "D:/samples/richness/wcvp/parquet"

sample_data <- open_dataset(
  sources = path, 
  format = "parquet"
)

results <- sample_data %>% 
  filter(sp <= 60000) %>% 
  filter(continent == "GLOBAL") %>% 
  collect() %>% 
  setDT() 

permutation_test100 <- function(results, n ) {
  p_frame <- data.frame(matrix(nrow = 1, ncol = 2))
  names(p_frame) <- c("wilcox_p", "sp")
  seq <- seq(n, n + 100,1)
  p_data <- results %>% 
    ungroup() %>% 
    filter(sp %in% seq) %>% 
    dplyr::select(cor.sp_gwr_mean_gwr_mean, sp)  %>% 
    mutate(sp = ifelse(!sp == n, n+100, sp)) 
  ttest <- perm.wilcox.test(cor.sp_gwr_mean_gwr_mean ~ sp, p_data, R = 999, alternative = "less")
  wtest <- wilcox.test(cor.sp_gwr_mean_gwr_mean ~ sp, p_data,alternative = "less" )
  p_frame$perm_p <- ttest$perm.p.value
  p_frame$wilcox_p <- wtest$p.value
  
  p_frame$g1 <-ttest$estimate[1]
  
  p_frame$g1 <-ttest$estimate[1]
  p_frame$g2 <-ttest$estimate[2]
  
  p_frame$sp <- paste(unique(p_data$sp)[1],"-",unique(p_data$sp[2]))
  return(p_frame)
}

permutation_test <- function(results, n ) {
  p_frame <- data.frame(matrix(nrow = 1, ncol = 2))
  names(p_frame) <- c("p", "sp")
  p_data <- results2 %>% 
    filter(sp == n | sp == n + 100) %>% 
    #mutate(sp = ifelse(!sp == n, n+100, sp)) %>% 
    dplyr::select(cor.sp_gwr_mean, sp) 
  ttest <- perm.wilcox.test(cor.sp_gwr_mean ~ sp, data = p_data,  R = 999, alternative = "less")
  p_frame$perm_p <- ttest$perm.p.value
  p_frame$p <-ttest$p.value
  p_frame$sp <- paste(unique(p_data$sp)[1],"-",unique(p_data$sp[2]))
  return(p_frame$perm_p)
}

#as.data.table(results2)[, perm.wilcox.test(cor.sp_gwr_mean ~ sp, data = p_data,  R = 999, alternative = "less"), by = c("continent", "sp")]

permutation_test100(results, n = 605)
results$pseudo.r.squared
n <- 500
permutation_test_cont <- function(results, n ) {

  p_data <- results %>% 
    group_by(continent) %>% 
    filter(sp == n | sp == n + 100) %>% # sp == ceiling(n + n*0.25)
    dplyr::select(pseudo.r.squared, sp, continent) 
  
  unique_sp <- unique(p_data$sp)
  
  test_data <- p_data %>%
    mutate(steps = paste(unique_sp[1], "-", unique_sp[2]),
           sp = as.factor(sp)) %>%
    group_by(continent) %>%
    reframe(perm_p = coin::pvalue(oneway_test(pseudo.r.squared ~ sp,  
                                        distribution = approximate(nresample = 999),
                                                                   #parallel = "snow",
                                                                   #cl = cl), 
                                        alternative = "less"))[1], 
            wilcox_p = wilcox.test(pseudo.r.squared ~ sp,  alternative = "less")$p.value,
            n = n,
            steps = unique(steps), 
            continent = unique(continent)) 
  
  return(test_data)
}

t_test_cont <- function(results, n ) {
  
  p_data <- results %>% 
    group_by(continent) %>% 
    filter(sp == n | sp == n + 100) %>% 
    #mutate(sp = ifelse(!sp == n, n+100, sp)) %>% 
    # sp == ceiling(n + n*0.25)
    dplyr::select(cor.sp_gwr_mean, sp, continent) 
  
  unique_sp <- unique(p_data$sp)
  
  test_data <- p_data %>%
    mutate(steps = paste(unique_sp[1], "-", unique_sp[2]),
           sp = as.factor(sp)) %>%
    group_by(continent) %>%
    reframe(wilcox_p = wilcox.test(cor.sp_gwr_mean ~ sp, alternative = "greater")$p.value,
            n = n,
            steps = unique(steps), 
            continent = unique(continent)) 
  return(test_data)
}

plan(multisession, workers = 6)

permutation_test_cont(results = results, n = 100)
permutation_test100(results = results, n = 501)


start <- Sys.time()
a <- future_lapply(seq(15,5000, 1), FUN = permutation_test_cont, results = results, future.seed = T)
#b <- lapply(1:1500, FUN = permutation_test100, results = results)
Sys.time() - start

xx <- do.call(rbind,a)


plot(
  xx$wilcox_p ~ xx$n)
plot(xx$perm_p ~   xx$wilcox_p )
lines(x=c(0,5000), y= c(0.05, 0.05), lwd =3, col = "red")
boxplot(xx$wilcox_p ~ xx$continent)
hist(xx$wilcox_p)
hist(xx$perm_p)


test <- xx %>% 
  group_by(continent) %>% 
  filter(wilcox_p >= 0.05) %>% 
  mutate(seq = case_when(lag(n) == n-1 ~ "yes", 
                         TRUE ~ "no")) %>% 
  mutate(consecutive = sequence(rle(as.character(seq))$lengths))

aa <- test %>% 
  filter(seq == "yes") %>% 
  slice_max(consecutive, n = 100) 


test %>% 
  filter(seq == "yes") %>% 
  filter(consecutive >= 10) %>% 
  slice_min(n, n = 1)

aa <- unique(test$n)
plot(aa)
boxplot(aa)

test2 <- results %>% 
  filter(sp %in% c(test$n, (test$n * 0.25))) %>% 
  group_by(sp, continent) %>% 
  reframe(mean = mean(cor.sp_gwr_mean, na.rm = T))
