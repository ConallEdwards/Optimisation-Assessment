#Library for finance analysis
library(tidyquant) 

#Library for handling portfolio optimisation
library(timetk)

#Libraries for data transformation
library(dplyr)
library(tidyr)

#Libraries for optimisation
library(GA)
library(mco)
library(metaheuristicOpt)
library(DEoptim)
library(pso)

#Libraries for visualisation
library(ggplot2)

#Library for statistical testing
library(stats)

#Setting seed for reproducability
set.seed(100)

load_data <- function(){
  #Loading in top 10 stocks
  ABBV_stock <- read.csv("ABBV.csv", header=TRUE)
  ABC_stock <- read.csv("ACB.csv", header=TRUE)
  BUD_stock <- read.csv("BUD.csv", header=TRUE)
  CGC_stock <- read.csv("CGC.csv", header=TRUE)
  CRON_stock <- read.csv("CRON.csv", header=TRUE)
  MJ_stock <- read.csv("MJ.csv", header=TRUE)
  MO_stock <- read.csv("MO.csv", header=TRUE)
  TAP_stock <- read.csv("TAP.csv", header=TRUE)
  TLRY_stock <- read.csv("TLRY.csv", header=TRUE)
  WEEDTO_stock <- read.csv("WEED.TO.csv", header=TRUE)
  
  #Add a symbol column to state which stock is which 
  ABBV_stock$Symbol <- "ABBV"
  ABC_stock$Symbol <- "ABC"
  BUD_stock$Symbol <- "BUD"
  CGC_stock$Symbol <- "CGC"
  CRON_stock$Symbol <- "CRON"
  MJ_stock$Symbol <- "MJ"
  MO_stock$Symbol <- "MO"
  TAP_stock$Symbol <- "TAP"
  TLRY_stock$Symbol <- "TLRY"
  WEEDTO_stock$Symbol <- "WEED.TO"
  
  #Binding all stock into one dataset
  stock_prices <- rbind(ABBV_stock, ABC_stock, BUD_stock, CGC_stock, CRON_stock, MJ_stock, MO_stock, TAP_stock, TLRY_stock, WEEDTO_stock)
  
  #Converting dates to proper type and format
  dates <- stock_prices$Date
  dates <- data.frame(as.Date(dates, format = "%Y-%m-%d"))
  stock_prices["Date"] <- dates
  
  
  #Selecting only the date, closing price, and the name of the stock 
  stock_prices_long <- stock_prices[c("Date", "Close", "Symbol")]
  stock_prices_long
  
  
  #Converting to wide format
  stock_prices_wide <- stock_prices_long %>%
    spread(Symbol, value = Close)
  
  stock_prices_wide <- na.omit(stock_prices_wide)
  
  stock_prices_wide <- subset(stock_prices_wide, select = -c(Date))
  
  return (stock_prices_wide)
  
}



#Function for optimisation methods that maximise (GA)
getFitnessGA <- function(weights){
  
  #Standardising the 
  weights <- weights/sum(weights)
  
  #Calculating portfolio returns 
  portfolio_return <- sum(weights * colMeans(stock_data))
  
  #Calculating portfolio variance
  portfolio_risk <- sqrt(t(weights) %*% cov(stock_data) %*% weights)
  
  #Calculating Sharpe ratio
  sharpeRatio <- portfolio_return / portfolio_risk
  
  
  return(sharpeRatio)
}




#Function for optimisation methods that minimise (PSO, DE)
getFitness <- function(weights){
  
  weights <- weights/sum(weights)
  
  #Calculating portfolio returns 
  portfolio_return <- sum(weights * colMeans(stock_data))
  
  #Calculating portfolio variance
  portfolio_risk <- sqrt(t(weights) %*% cov(stock_data) %*% weights)
  
  #Calculating Sharpe ratio
  sharpeRatio <- portfolio_return / portfolio_risk
  
  
  return(-sharpeRatio)
}



runGA <- function(){
  results <- ga(type ="real-valued", fitness = getFitnessGA, lower = rep(0, 10), upper = rep(1, 10), 
                monitor = TRUE, maxiter = 100, keepBest = TRUE)
  return(results)
}


runPSO <- function(){
  results <- psoptim(rep(NA, 10), fn = getFitness, lower=rep(0, 10), upper = rep(1, 10), 
                     control = c(maxit = 100, REPORT = 1, trace = 1, trace.stats=1))
  return(results)
}


runDE <- function(){
  results <- DEoptim(fn = getFitness, lower = rep(0, 10), upper = rep(1, 10),  
                     DEoptim.control(itermax = 100))
  return(results)
}


#Loading dataset
stock_data <- load_data()

#Getting benchmark value
benchmark_solution <- c(0.1, 0.1, 0.1, 0.1, 0.1,0.1, 0.1, 0.1, 0.1, 0.1)
benchmark_sharpe_ratio = getFitnessGA(benchmark_solution)


#Running experiments
GAresults <- runGA()
PSOresults <- runPSO()
DEresults <- runDE()


#Getting the optimum solutions
GAOptim <- GAresults@solution / sum(GAresults@solution)
GAOptim


PSOOptim <- PSOresults$par / sum(PSOresults$par)
PSOOptim

DEOptim <- DEresults$optim$bestmem / sum(DEresults$optim$bestmem)
DEOptim


#Calculate fitness values for GA and PSO
GA_fitnesses <- c()
PSO_fitnesses <- c()

for (i in c(1:100)){
  GA_fit_value <- getFitnessGA(GAresults@bestSol[[i]][1,]) #change this to calculate sharpe ratio here! 
  GA_fitnesses <- append(GA_fitnesses, GA_fit_value)
  
  PSO_fit_value <- min(PSOresults$stats$f[[i]][1])
  PSO_fitnesses  <- append(PSO_fitnesses, PSO_fit_value)
  
}


#Calculating standard deviations of each fitness scores
sd_ga_value <- sd(GA_fitnesses)
sd_pso_value <- sd(PSO_fitnesses)
sd_de_value <- sd(DEresults$member$bestvalit)


#Converting arrays to dataframes 
GA_fitnesses <- as.data.frame(GA_fitnesses)
PSO_fitnesses <- as.data.frame(PSO_fitnesses)
DE_fitnesses <- as.data.frame(DEresults$member$bestvalit)

#Concatenating all fitnesses into one dataframe
all_fitnesses <- cbind(GA_fitnesses, PSO_fitnesses, DE_fitnesses)
colnames(all_fitnesses) <- c("GA_Fitness", "DE_Fitness", "PSO_Fitness")


#Plotting all fitnesses in one graph
ggplot(all_fitnesses, aes(x = c(1:100))) +
  geom_line(aes(y = GA_Fitness, color = "red")) +
  geom_line(aes(y = -DE_Fitness, color = "blue")) +
  geom_line(aes(y = -PSO_Fitness, color = "darkgreen")) +
  geom_errorbar(aes(ymin = GA_Fitness - sd_ga_value, ymax = GA_Fitness + sd_ga_value, color = "red")) +
  geom_errorbar(aes(ymin = -DE_Fitness - sd_de_value , ymax = -DE_Fitness + sd_de_value, color = "blue")) +
  geom_errorbar(aes(ymin = -PSO_Fitness - sd_pso_value, ymax = -PSO_Fitness + sd_pso_value, color = "darkgreen")) +
  xlab("Number of generations") + 
  ylab("Fitness value") +
  ggtitle("Best fitness value achieved by each optimisation technique over 100 generations") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.85, 0.8)) +
  labs(colour = "Optimisation Algorithm") + 
  scale_color_manual(labels = c("DE", "GA", "PSO"), values = c("blue", "red", "darkgreen"))


#Plotting the GA
plot(x = c(1:GAresults@iter) , 
     y = GAresults@summary[,"mean"],
     main = "Showing how the GA evolves to get the best Sharpe ratio over 100 generations",
     xlab = "Number of generations",
     ylab = "Fitness value",
     type = "l")


#Plotting the PSO
plot(x = c(1:100), 
     y = -PSO_fitnesses[,1],
     main = "Showing how the PSO evolves to get the best Sharpe ratio over 100 generations",
     xlab = "Number of generations",
     ylab = "Fitness value",
     type = "l"
     )


#Plotting the DE
plot(x = c(1:100), 
     y = -DE_fitnesses[,1],
     main = "Showing how the DE evolves to get the best Sharpe ratio over 100 generations",
     xlab = "Number of generations",
     ylab = "Fitness value",
     type = "l"
)



#Apply t-tests to fitness scores achieved
t.test(x = GA_fitnesses,mu = benchmark_sharpe_ratio, alternative = "greater")
t.test(x = PSO_fitnesses, mu = benchmark_sharpe_ratio, alternative = "greater")
t.test(x = DE_fitnesses, mu = benchmark_sharpe_ratio, alternative = "greater")
