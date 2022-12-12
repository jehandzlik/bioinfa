library(ggplot2)
library(plotly)
library(ggfortify)

output.dir <- "/Users/leukosnanos/Desktop/Bioinfromatics_for_Computer_Scientistis/Plots/SVD"
setwd(output.dir)

## SVD

calcEigen <- function(df, type){
  
  cov.type <- cov(df[which(df$type == type),1:3])
  eigen.vectors.values <- eigen(cov.type)
  eigen.vec <- eigen.vectors.values$vectors
  eigen.val <- eigen.vectors.values$values
  return(list(eigen.vec, eigen.val))
  
}

createData <- function(n, centered.name){
  
  set.seed(7)
  ## Generate two white 3D data sets, one centered around
  ## zero and one not.
  x.centered = rnorm(n, mean = 0, sd = 1)
  y.centered = rnorm(n, mean = 0, sd = 1)
  z.centered = rnorm(n, mean = 0, sd = 1)

  df <- data.frame(x = x.centered,
                   y = y.centered,
                   z = z.centered, type = centered.name)
  
}

set.seed(7)
## Generate two white 3D data sets, one centered around
## zero and one not.
n = 1000
x.centered = rnorm(n, mean = 0, sd = 1)
y.centered = rnorm(n, mean = 0, sd = 1)
z.centered = rnorm(n, mean = 0, sd = 1)
x.not.centered = x.centered + 7
y.not.centered = y.centered + 12
z.not.centered = z.centered + 13

centered.name = "centered"
non.centered.name = "non.centered"
centered.transformed.name = "centered.transformed"
non.centered.transformed.name = "non.centered.transformed"

df <- data.frame(x = x.centered,
                 y = y.centered,
                 z = z.centered, type = centered.name)
df <- rbind(df, data.frame(x = x.not.centered,
                           y = y.not.centered,
                           z = z.not.centered,
                           type = non.centered.name))


# transformed <- matrix(c(4,0,0,1), nrow = 2, ncol = 2) %*% t(df[,1:2])
# plot(transformed[1,], transformed[2,])
# plot(df[,1], df[,2])
# 
# cov(t(transformed))

## Plot the centered and not-centered data

plot_ly(x=df$x, y=df$y, z=df$z, type="scatter3d", mode="markers", color=df$type, size = 0.1)

pca <- prcomp(df[which(df$type == centered.name | df$type == non.centered.name),1:3], 
              scale. = TRUE, center = TRUE)
autoplot(pca)

# ggplot(df, aes(x = x, y = y, color = type)) +
#   geom_point()

## The covariance matrix is the same for the
## centered and non centered data.
cov.centered <- cov(df[which(df$type == centered.name),1:3])
cov.non.centered <- cov(df[which(df$type == non.centered.name),1:3])

## If you do eigen on df directly, you get the error
## "non-square matrix in 'eigen'"
results <- calcEigen(df, centered.name)
eigen.vec <- results[1][[1]]
eigen.val <- as.vector(results[2][[1]])

# p1 <- ggplot(df, aes(x = x, y = y, group = type, color = type)) +
#   geom_point() +
#   geom_segment(x = 0, y = 0, xend = eigen.vec[1,1]*eigen.val[1],
#                yend = eigen.vec[2,1]*eigen.val[1], col = "blue", size = 3) +
#   geom_segment(x = 0, y = 0, xend = eigen.vec[1,2]*eigen.val[2],
#                yend = eigen.vec[2,2]*eigen.val[2], col = "blue", size = 3)

#ggsave(filename = "data.pdf", p1, width = 16, height = 10, units = "in")


############################################
## Apply transformation onto the white data.
############################################

## Create a transformation matrix
#T = matrix(c(1,2,0.5,1,2,2,0.4,0.1,1.2), nrow = 3, ncol = 3) ## very strange transformation
T = matrix(c(1,2,0.5,4,2,2,0.4,0.1,1.2), nrow = 3, ncol = 3)

#T = matrix(c(2,3,10,1,1,1,1,1,1), nrow = 3, ncol = 3)

#T = matrix(c(2,2,2,3,3,3,1,1,1), nrow = 3, ncol = 3)

## Transform the data
centered.transformed <- T %*% t(as.matrix(df[which(df$type == centered.name),1:3]))
non.centered.transformed <- T %*% t(as.matrix(df[which(df$type == non.centered.name),1:3]))


T = matrix(c(3,0,2,2), nrow = 2, ncol = 2)

eigen(T)

twoD.transformed <- T %*% t(as.matrix(df[which(df$type == centered.name),1:2]))


cov(t(twoD.transformed))

eigen()


plot(twoD.transformed[1,], twoD.transformed[2,])

cov(t(twoD.transformed))

Tsym = matrix(c(3,1,1,2), nrow = 2, ncol = 2)
twoD.transformed.sym <- Tsym %*% t(as.matrix(df[which(df$type == centered.name),1:2]))

plot(twoD.transformed.sym[1,], twoD.transformed.sym[2,])
cov(t(twoD.transformed.sym))

## The covariance matrix of transformed white data is
## equal to the transformation matrix multiplied by itself, T*T'
cov(t(centered.transformed))
T %*% t(T)

cov(t(non.centered.transformed))

## Add transformed data to data
df <- rbind(df, data.frame(x = centered.transformed[1,],
                           y = centered.transformed[2,],
                           z = centered.transformed[3,],
                           type = centered.transformed.name))

df <- rbind(df, data.frame(x = non.centered.transformed[1,],
                           y = non.centered.transformed[2,],
                           z = non.centered.transformed[3,],
                           type = non.centered.transformed.name))

plot_ly(x=df$x, y=df$y, z=df$z, type="scatter3d", mode="markers", color=df$type, size = 0.1)

df.prt <- df[which(df$type == "centered" | df$type == "centered.transformed"),]

plot_ly(x=df.prt$x, y=df.prt$y, z=df.prt$z, type="scatter3d", mode="markers", color=df.prt$type, size = 0.1)

cov.transformed <- cov(t(non.centered.transformed)) %*% t(as.matrix(df[which(df$type == centered.name),1:3]))

cov(t(cov.transformed))

df.prt <- rbind(df.prt, data.frame(x = cov.transformed[1,],
                           y = cov.transformed[2,],
                           z = cov.transformed[3,],
                           type = "cov.transformed"))

plot_ly(x=df.prt$x, y=df.prt$y, z=df.prt$z, type="scatter3d", mode="markers", color=df.prt$type, size = 0.1)


# results <- calcEigen(df, centered.transformed.name)
# eigen.vec <- results[1][[1]]
# eigen.val <- as.vector(results[2][[1]])

#######################
## Perform PCA manually
#######################

df.prt <- df[which(df$type == "centered.transformed"),]

cov.type <- cov(df.prt[,1:3])
eigen.vectors.values <- eigen(cov.type)
eigen.vec <- eigen.vectors.values$vectors
eigen.val <- eigen.vectors.values$values


eigen.vectors.values <- eigen(T)
eigen.vec.T <- eigen.vectors.values$vectors
eigen.val.T <- eigen.vectors.values$values

#inv <- solve(eigen.vec)
pca.manual <- t(as.matrix(eigen.vec[,1:2])) %*% t(as.matrix(df.prt[,1:3])) #project the data into the 2D space

pca.manual <- t(as.matrix(eigen.vec.T[,1:2])) %*% t(as.matrix(df.prt[,1:3])) #project the data into the 2D space
#pca.manual <- t(as.matrix(inv[,1:2])) %*% t(as.matrix(df.prt[,1:3])) #project the data into the 2D space

pca.manual <- t(pca.manual)

scaling.manual <- c(sd(pca.manual[,1]) * sqrt(nrow(df.prt)), sd(pca.manual[,2]) * sqrt(nrow(df.prt)))
plot(pca.manual[,1]/scaling.manual[1], pca.manual[,2]/scaling.manual[2])

## Automatic PCA plot
pca <- prcomp(df.prt[,1:3], 
              scale. = TRUE, center = TRUE)
autoplot(pca)

plot(pca.manual[,1], pca.manual[,2])


df.prt <- rbind(df.prt, data.frame(x = pca.manual[,1], y = pca.manual[,2], z = 0, type = "reduced.T"))
prev <- df.prt[which(df.prt$type == "reduced"),]
df.prt <- rbind(df.prt, prev)

df.prt <- df.prt[which(df.prt$type == "reduced" | df.prt$type == "reduced.T"),]


plot_ly(x=df.prt$x, y=df.prt$y, z=df.prt$z, type="scatter3d", mode="markers", color=df.prt$type, size = 0.1)

## Plot eigenvectors
m = 50
PC1 <- eigen.vec[,1] * m
PC2 <- eigen.vec[,2] * m
PC3 <- eigen.vec[,3] * m

PC1.T <- eigen.vec.T * m
PC2.T <- eigen.vec.T[,2] * m
PC3.T <- eigen.vec.T[,3] * m

k = 50
## Plot all the points and the 2D surface that the two eigenvectors form
plot_ly(df.prt, color = df.prt$type, size = 0.01) %>%
  add_markers(x = ~x, y = ~y, z = ~z) %>%
  add_mesh(x = c(0, PC1[1], PC2[1]),
           y = c(0, PC1[2], PC2[2]),
           z = c(0, PC1[3], PC2[3]),
           opacity = 0.3) %>%
  add_mesh(x = c(0, PC1.T[1], PC2.T[1]),
           y = c(0, PC1.T[2], PC2.T[2]),
           z = c(0, PC1.T[3], PC2.T[3]),
           opacity = 0.3)
  # %>%
  # add_mesh(x = c(0, eigen.vec[1,1], eigen.vec[1,2])*m,
  #           y = c(0, eigen.vec[1,2], eigen.vec[2,2]*m),
  #           z = c(0, eigen.vec[1,3], eigen.vec[3,2]*m))

colnames(data) <- c("a", "b", "c", "color")

## Find a line equation for the the first eigenvector
lmn = eigen.vec[,1] - c(0,0,0)
l <- lmn[1]
m <- lmn[2]
n <- lmn[3]

bcd = eigen.vec[,2] - c(0,0,0)
b <- bcd[1]
c <- bcd[2]
d <- bcd[3]

efg = eigen.vec[,3] - c(0,0,0)
e <- efg[1]
f <- efg[2]
g <- efg[3]

n = 1000

a <- runif(n)
# y <- runif(n)
# z <- runif(n)


data <- data.frame()
data <- rbind(data, data.frame(x = (a-lmn[1])/l, y = (a-lmn[2])/m, z = (a-lmn[3])/n, color = rep(1,n)))
data <- rbind(data, data.frame(x = (a-bcd[1])/b, y = (a-bcd[2])/c, z = (a-bcd[3])/d, color = rep(2,n)))
data <- rbind(data, data.frame(x = (a-efg[1])/e, y = (a-efg[2])/f, z = (a-efg[3])/g, color = rep(3,n)))
#plot_ly()
#data <- read.csv('https://raw.githubusercontent.com/plotly/datasets/master/3d-line1.csv')

data$color <- as.factor(data$color)


plot_ly(df.prt, color = df.prt$type, size = 0.01) %>%
  add_markers(x = ~x, y = ~y, z = ~z) %>%
  add_mesh(x = c(0, PC1[1], PC2[1]), 
           y = c(0, PC1[2], PC2[2]),
           z = c(0, PC1[3], PC2[3]), 
           opacity = 0.3) %>%
  add_lines(x = data[,1], 
           y = data[,2],
           z = data[,3])


plot_ly(data, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines',
        
        opacity = 1, line = list(width = 6, color = ~color, reverscale = FALSE))

plot_ly(data, color = data$color) %>%
  add_lines(x = ~x, 
            y = ~y,
            z = ~z)

cov.type <- cov(df[,1:3])
eigen.vectors.values <- eigen(cov.type)
eigen.vec <- eigen.vectors.values$vectors
eigen.val <- eigen.vectors.values$values

pca.manual <- t(as.matrix(eigen.vec[,1:2])) %*% t(as.matrix(df[,1:3])) #project the data into the 2D space

pca.manual <- t(pca.manual)

scaling.manual <- c(sd(pca.manual[,1]) * sqrt(nrow(df)), sd(pca.manual[,2]) * sqrt(nrow(df)))
plot(pca.manual[,1]/scaling.manual[1], pca.manual[,2]/scaling.manual[2])

plot(pca.manual[,1], pca.manual[,2])

data
svd(data)

## Project the data into the new axis (2d surface formed by the two eigenvectors)
pca.manual <- t(as.matrix(pca.rotation[,1:2])) %*% t(as.matrix(data))

pca.manual <- t(pca.manual)

scaling.manual <- c(sd(pca.manual[,1]) * sqrt(nrow(data)), sd(pca.manual[,2]) * sqrt(nrow(data)))
plot(pca.manual[,1]/scaling.manual[1], pca.manual[,2]/scaling.manual[2])
## OR
autoplot(pca.prcomp)


scaling <- pca.prcomp$sdev[1:2] * sqrt(nrow(data))

s.eigen <- eigen(cov(data))

pc1 <- rowSums(t(t(sweep(data, 2 ,colMeans(data))) * s.eigen$vectors[,1] * -1) / scaling[1])
pc2 <- rowSums(t(t(sweep(data, 2, colMeans(data))) * s.eigen$vectors[,2]) / scaling[2])

#Collect the PCs in a data.frame and plot using ggplot (loaded when ggfortify was loaded).

df <- data.frame(pc1, pc2, c(rep('Apprentice', 20), rep('Pilot', 20)))
colnames(df) <- c('PC1', 'PC2', 'Group')

ggplot(df, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point()

############################################
## Now work only on centered and scaled data
############################################

data <- df[which(df$type == centered.transformed.name),1:3]

#cov(data)
## Perform the R PCA
pca.prcomp <- prcomp(data)
pca.prcomp.rotation <- pca.prcomp$rotation
summary(pca.prcomp)

m = 50
PC1 <- pca.prcomp.rotation[,1] * -m
PC2 <- pca.prcomp.rotation[,2] * -m
#PC3 <- pca.prcomp.rotation[,3] * -m

## Plot all the points and the 2D surface that the two eigenvectors form
plot_ly(df, color = df$type, size = 0.01) %>%
  add_markers(x = ~x, y = ~y, z = ~z) %>%
  add_mesh(x = c(0, PC1[1], PC2[1]), 
           y = c(0, PC1[2], PC2[2]),
           z = c(0, PC1[3], PC2[3]), 
           opacity = 0.3)


#######################
## Perform PCA manually
#######################

data
svd(data)

## Project the data into the new axis (2d surface formed by the two eigenvectors)
pca.manual <- t(as.matrix(pca.rotation[,1:2])) %*% t(as.matrix(data))

pca.manual <- t(pca.manual)

scaling.manual <- c(sd(pca.manual[,1]) * sqrt(nrow(data)), sd(pca.manual[,2]) * sqrt(nrow(data)))
plot(pca.manual[,1]/scaling.manual[1], pca.manual[,2]/scaling.manual[2])
## OR
autoplot(pca.prcomp)


scaling <- pca.prcomp$sdev[1:2] * sqrt(nrow(data))

s.eigen <- eigen(cov(data))

pc1 <- rowSums(t(t(sweep(data, 2 ,colMeans(data))) * s.eigen$vectors[,1] * -1) / scaling[1])
pc2 <- rowSums(t(t(sweep(data, 2, colMeans(data))) * s.eigen$vectors[,2]) / scaling[2])

#Collect the PCs in a data.frame and plot using ggplot (loaded when ggfortify was loaded).

df <- data.frame(pc1, pc2, c(rep('Apprentice', 20), rep('Pilot', 20)))
colnames(df) <- c('PC1', 'PC2', 'Group')

ggplot(df, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point()

##########################################
## Combine centered and not centered data
##########################################

pca <- prcomp(df[which(df$type == centered.name | df$type == non.centered.transformed.name),1:3], 
              scale. = TRUE, center = TRUE)
autoplot(pca)

autoplot(pca, scale = TRUE, center = TRUE)

# df.test <- data.frame(x = rnorm(10), y = rnorm(10), z = rnorm(10))
# planDf <- data.frame(x = rep(range(df.test$x), 2), y = rep(range(df.test$y), each = 2), z = mean(df.test$z))
# 
# plot_ly(df.test) %>%
#   add_markers(x = ~x, y = ~y, z = ~z) %>%
#   add_mesh(x = ~x, y = ~y, z = ~z, data = planDf, opacity = 0.3)

# # initiate a line shape object
#
# line <- list(
#   type = "line",
#   line = list(color = "pink"),
#   xref = "x",
#   yref = "y",
#   zref = "z"
# )
#
#
# lines <- list()
#
# line[["x0"]] <- 0
# line[["y0"]] <- 0
# line[["z0"]] <- 0
# line[c("x1", "y1", "z1")] <- c(eigen.vec[1,1], eigen.vec[2,1], eigen.vec[3,1])
# lines <- c(lines, list(line))
#
# # for (i in c(0, 3, 5, 7, 9, 13)) {
# #   line[["x0"]] <- i
# #   line[["x1"]] <- i + 2
# #   line[c("y0", "y1")] <- sin(i + 1)
# #   lines <- c(lines, list(line))
# # }
#
# fig <- layout(fig, title = 'Highlighting with Lines', shapes = lines)
#
# fig.add_trace(go.Scatter(x=eigen.vec[1,1], y=eigen.vec[2,1], z = eigen.vec[3,1],
#                          mode='lines', marker=dict(symbol="line-ne"),
#                          line=dict(color='red',width=2), name='week'))
#
#
#
# ## Plot the data and the transformed data
# p2 <- ggplot(df, aes(x = x, y = y, group = type, color = type)) +
#   geom_point(size = 0.3)+
#   geom_segment(x = 0, y = 0, xend = eigen.vec[1,1]*eigen.val[1],
#                yend = eigen.vec[2,1]*eigen.val[1], col = "blue", size = 2) +
#   geom_segment(x = 0, y = 0, xend = eigen.vec[1,2]*eigen.val[2],
#                yend = eigen.vec[2,2]*eigen.val[2], col = "blue", size = 2)
# #geom_segment(x = 0, y = 0, xend = T[1,1], yend = T[2,1], col = "red") +
# #geom_segment(x = 0, y = 0, xend = T[2,1], yend = T[2,2], col = "red")
# #geom_segment(x = 0, y = 0, xend = 17, yend = 10, col = "red") +
# #geom_segment(x = 0, y = 0, xend = 10, yend = 8, col = "red")
#
# ggsave(filename = "transformed.pdf", p2, width = 16, height = 10, units = "in")
# #
# cov(centered.transformed)
# T %*% t(T)


# #######################
# ## PCA
# #######################

# data = centered.transformed.name
# #data = non.centered.transformed.name
# #data = centered.name

# data <- df[which(df$type == data),1:3]
# 
# plot(data[,1], data[,2])
# plot(data[,2], data[,3])
# plot(data[,1], data[,3])
# 
# cov(data)
# 
# pca <- prcomp(data, scale. = T, center = T)
# pca <- prcomp(data)
# pca$sdev^2/sum(pca$sdev^2)
# autoplot(pca)




# #new.coordinates <- scale(solve(pca$rotation)) %*% as.matrix(t(df.part))
#
# new.coordinates <- eigen.vec %*% as.matrix(t(centered.transformed.df))
#
# new.coordinates <- eigen.vec %*% new.coordinates
#
# dot.product <- eigen.vec[,1] %*% as.matrix(t(centered.transformed.df))
#
# new.coordinates.PC1 <- eigen.vec[,1] * dot.product
#
# ## project new data onto the PCA space
# ## new.coordinates <- scale(df.part, pca$center, pca$scale) %*% pca$rotation
#
# df.backup <- df
# df <- df.backup
#
# df <- rbind(df, data.frame(x = new.coordinates[1,],
#                            y = new.coordinates[2,],
#                            type = "new.coordinates"))
#
# ggplot(df, aes(x = x, y = y, group = type, color = type)) +
#   geom_point(size = 0.9) +
#   geom_segment(x = 0, y = 0, xend = eigen.vec[1,1]*eigen.val[1],
#                yend = eigen.vec[2,1]*eigen.val[1], col = "blue", size = 1) +
#   geom_segment(x = 0, y = 0, xend = eigen.vec[1,2]*eigen.val[2],
#                yend = eigen.vec[2,2]*eigen.val[2], col = "blue", size = 1)
#

n = 100
x.centered = rnorm(n, mean = 0, sd = 1)

a <- seq(1,10,by = 1)
sd(a)
b <- seq(1,20,by = 2)
sd(b)



############################################
## Eigenvectors of a non-symmetric matrix
## pointing towards the highes variance
############################################

centered.name = "centered"
non.centered.name = "non.centered"
centered.transformed.name = "centered.transformed"
non.centered.transformed.name = "non.centered.transformed"

df <- createData(1000, centered.name)

data <- as.matrix(df[which(df$type == centered.name),1:2])

## Symmetric transformation matrix
A = matrix(c(3,2,2,1), nrow = 2, ncol = 2)
## OR
## Non-symmetric transformation matrix
A = matrix(c(5,1,2,2), nrow = 2, ncol = 2)

transformed.data <- t(A %*% t(data[,1:2]))

eigen.vectors.values <- eigen(A)
eigen.vec.A <- eigen.vectors.values$vectors
eigen.val.A <- eigen.vectors.values$values

eigen.vectors.values <- eigen(A %*% t(A))
eigen.vec.AA <- eigen.vectors.values$vectors
eigen.val.AA <- eigen.vectors.values$values

data.df <- data.frame(x = data[,1], y = data[,2], type = "centered")
data.df <- rbind(data.df, data.frame(x = transformed.data[,1], y = transformed.data[,2], type = "transformed"))

ggplot(data.df, aes(x = x, y = y, color = type)) +
  geom_point() +
  geom_abline(intercept = 0, slope = eigen.vec.A[2,1]/eigen.vec.A[1,1]) +
  geom_abline(intercept = 0, slope = eigen.vec.AA[2,1]/eigen.vec.AA[1,1]) +
  ylim(c(-10,10))

## Create a transformation matrix

A = matrix(c(3,0,2,2), nrow = 2, ncol = 2)

eigen.vectors.values <- eigen(A %*% t(A))
eigen.vec.AA <- eigen.vectors.values$vectors
eigen.val.AA <- eigen.vectors.values$values


A %*% eigen.vec.A


eigen.vectors.values <- eigen(A)
eigen.vec.A <- eigen.vectors.values$vectors
eigen.val.A <- eigen.vectors.values$values

transformed.data <- t(A %*% t(data[,1:2]))

plot(transformed.data[,1], transformed.data[,2], xlim = c(-10,10), ylim = c(-10,10))
abline(a = 0, b = eigen.vec.A[2,1]/eigen.vec.A[1,1])
abline(a = 0, b = eigen.vec.A[2,2]/eigen.vec.A[1,2])

## Calculate the eignevectors of the transformed data
## or the direction of the maximum variance in the data
## A) using covariance matrix


plot(transformed.data[,1], transformed.data[,2], xlim = c(-10,10), ylim = c(-10,10))
abline(a = 0, b = eigen.vec.AA[2,1]/eigen.vec.AA[1,1])
abline(a = 0, b = eigen.vec.AA[2,2]/eigen.vec.AA[1,2])

## B) using Av1 and Av2

transformed.data <- t(A %*% (A) %*% t(data[,1:2]))

plot(transformed.data[,1], transformed.data[,2], xlim = c(-10,10), ylim = c(-10,10))
abline(a = 0, b = eigen.vec.AA[2,1]/eigen.vec.AA[1,1])
abline(a = 0, b = eigen.vec.AA[2,2]/eigen.vec.AA[1,2])

twoD.transformed <- T %*% t()

eigen(cov(t(twoD.transformed)))

eigen.vec %*% T
