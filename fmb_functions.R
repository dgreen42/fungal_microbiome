# David C. Green
# 2023-03-28
# Microbiome functions
# ---------------------------------------------------------------------
# Libraries
library(tidyverse)
library(ggthemes)
library(patchwork)
library(igraph)
library(plotly)
library(rethinking)
library(beepr)

# ---------------------------------------------------------------------
# Sources

# ---------------------------------------------------------------------

predTest <- function(a, b, data) {
  x <- NULL
  for(i in 1:70) {
    n <- i + 1
    x[i] <- data[a,n]
  }
  y <- NULL
  for(i in 1:70) {
    n <- i + 1
    y[i] <- data[b,n]
  }
  df <- data.frame(x,y)
  linM <- lm(y ~ x, data = df)
  s <- summary(linM)
  r <- s$r.squared
  r <- format(r, digits = 2)
  plot(x,y,main = c(a, 'and', b), sub = r)
  abline(linM)
  return(r)
}


QUADfit <- function(a, b, data) {
  x <- NULL
  for(i in 1:70) {
    n <- i + 1
    x[i] <- data[a,n]
  }
  y <- NULL
  for(i in 1:70) {
    n <- i + 1
    y[i] <- data[b,n]
  }
  df <- data.frame(x,y)
  model <- quap(
    alist(
      y ~ dnorm(mu, sigma),
      mu <- a + b*x,
      a ~ dnorm(0, 0.2),
      b ~ dnorm(0, 0.5),
      sigma ~ dexp(1)),
    data = df
  )
  return(model)
}

correlationMatrix <- function(dataSet) {
  rMatrix <- matrix(nrow = nrow(dataSet), ncol = nrow(dataSet))
  for(i in 1:nrow(dataSet)) {
    for(j in 1:nrow(dataSet)) {
      x <- NULL
      for(l in 1:70) {
        n <- l + 1
        x[l] <- dataSet[i,n]
      }
      y <- NULL
      for(l in 1:70) {
        n <- l + 1
        y[l] <- dataSet[j,n]
      }
      df <- data.frame(x,y)
      linM <- lm(y ~ x, data = df)
      s <- summary(linM)
      r <- s$r.squared
      rMatrix[i,j] <- r
    }
  }
  beep('fanfare')
  return(rMatrix)
}

# For categorical variables across data frames
# x takes list of things to color by
# y takes unorganized list of things to assign colors to
# l takes output list containing x as the first column
# colors takes color gradient
vcolors <- function(x, y, colors) {
  l <- data.frame(x)
  collist <- colorRampPalette(colors)
  cols <- collist(nrow(x))
  for(i in 1: nrow(genL)) {
    tfl <- x[i,] == y
    for(j in 1:length(tfl)) {
      ifelse(tfl[j] == TRUE, l[j,2] <- cols[i], next)
    }
  }
  return(l[,2])
}

#colors for list of integers
vcolors2 <- function(x, colors) {
  collist <- colorRampPalette(colors)
  meanDF <- data.frame(x)
  meansS <- sort(x)
  colors <- collist(length(x))
  meanColDf <- data.frame(meansS, colors)
  for(i in 1:length(x)) {
    for(j in 1:length(x)) {
      ifelse(meanDF[i,1] == meanColDf[j,1], meanDF[i,2] <- meanColDf[j,2], next)
    }
  }
  return(meanDF$V2)
}

#average abundance
avgAbundance <- function(dataSet) {
  meanL <- NULL
  for(i in 1:nrow(dataSet)) {
    x <- NULL
    for(j in 1:length(dataSet)-1) {
      x[j] <- as.integer(dataSet[i,j+1])
    }
    meanL[i] <- mean(x)
  }
  return(meanL)
}