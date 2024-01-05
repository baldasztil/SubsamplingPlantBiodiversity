results <- read.csv("fullsamples_test.csv") %>% 
  mutate(sp = as.factor(sp))

install.packages("MKinfer")
library(MKinfer)

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

as.data.table(results2)[, perm.t.test(cor.sp ~ sp, data = p_data,  R = 999, alternative = "less"), by = c("id", "sp")]


results2 <- results %>% 
  filter(id %in% c("1", "2")) %>% 
  filter(sp == 100)

unique(results2$sp)

install.packages("coin")
library(coin)


permutation_test_cont <- function(results, n ) {

  p_data <- results %>% 
    group_by(id) %>% 
    filter(sp == n | sp == n + 100) %>% 
    dplyr::select(cor.sp, sp, id) 
  
  start <- Sys.time()
  test_data <- p_data %>% 
    group_by(id) %>% 
    reframe(perm_p = perm.t.test(cor.sp ~ sp, data = p_data,R = 999, 
                                 alternative = "less")$perm.p.value,
            steps = paste(unique(p_data$sp)[1],"-",unique(p_data$sp[2])),
            id = unique(id), 
            n = n)
  Sys.time() - start
  
  start <- Sys.time()
  test_data <- p_data %>%
    mutate(steps = paste(unique_sp[1], "-", unique_sp[2])) %>%
    group_by(id) %>%
    summarise(perm_p = perm.t.test(cor.sp ~ sp, data = ., R = 999, alternative = "less")$perm.p.value) %>%
    ungroup()
  Sys.time() - start
  
  start <- Sys.time()
  test_data <- p_data %>%
    mutate(steps = paste(unique_sp[1], "-", unique_sp[2]),
           sp = as.factor(sp)) %>%
    group_by(id) %>%
    reframe(perm_p = pvalue(oneway_test(cor.sp ~ sp, data = ., distribution = 
                                            approximate(nresample = 10000), 
                                   alternative = "less"))[1], 
             n = n,
            steps = unique(steps), 
            id = unique(id)) 
  Sys.time() - start
  
  return(test_data)
}


aa <- oneway_test(cor.sp ~ sp, data = p_data, distribution = approximate(nresample = 99999), 
            alternative = "less", conf.level = 0.95)
aa <- perm.t.test(cor.sp ~ sp, data = p_data, R = 999, alternative = "less")
pvalue(aa, level = 0.05)
pvalue_interval(aa)

aa@distribution@pvalue()
size(aa)

permutation_test_cont(results, n = 500)



results_sum <- results %>%  
  group_by(id) %>% 
  summarise()

permutation_test(results, 1500)


resultsxx <- results2 %>% 
  filter(sp <= 500) %>% 
  filter(id == "overall") %>% 
  mutate(sp = paste0("step",sp)) %>% 
  mutate(sp = as.factor(sp)) %>% 
  filter(!is.na(cor.sp)) %>% 
  dplyr::select(cor.sp, sp)

unique(results$sp)

a <- lapply(1:1500, FUN = permutation_test_cont, results = results)
#b <- lapply(1:1500, FUN = permutation_test100, results = results)

xx <- do.call(rbind,a)
xx$test <- "perm"
xx$id <- as.numeric(rownames(xx))

yy <- do.call(rbind,b)
yy$test <- "permall100"
yy$id <- as.numeric(rownames(yy))

plot(xx$perm_p ~ xx$id)
plot(yy$perm_p ~ yy$id)
boxplot(yy$perm_p)
boxplot(xx$perm_p)

zz <- rbind(xx, yy)

boxplot(zz$perm_p ~ zz$test)

test <- zz %>% 
  group_by(test) %>% 
  filter(test == "perm") %>% 
  filter(perm_p >= 0.05) %>% 
  slice_min(n = 20, order_by = id)

test2 <- results %>% 
  filter(sp %in% test$id) %>% 
  group_by(sp) %>% 
  summarise(mean = mean(cor.sp))
