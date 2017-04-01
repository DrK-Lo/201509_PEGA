#install.packages("ggraph")
#install.packages("tidyverse")
#install.packages("corrr")
#install.packages('dendextend')
library(ggraph)
library(igraph)
library(tidyverse)
library(corrr)
library("dendextend")
?graph_from_data_frame
head(highschool)

head(mtcars)
dim(mtcars)
mtcars %>% 
  correlate()

tidy_cors <- mtcars %>% 
  correlate() %>% 
  stretch()
head(tidy_cors)

d2<- dist(t(mtcars)) 
myFun <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
}
d3 <- myFun(d2)
# Next, we convert these values to an undirected graph object. 
# The graph is undirected because correlations do not have a 
# direction. For example, correlations do not assume cause or 
#effect. This is done using the igraph function, 
# graph_from_data_frame(directed = FALSE).

# Because, we typically don’t want to see ALL of the correlations, 
# we first filter() out any correlations with an absolute value 
# less than some threshold. For example, let’s include correlations 
# that are .3 or stronger (positive OR negative):

tidy_cors %>% filter(abs(r) > .3)

d3 %>% filter(abs(value) < 500)

graph_cors <- tidy_cors %>%
  filter(abs(r) > .3) %>%
  graph_from_data_frame(directed = FALSE)

graph_cors2 <- d3 %>%
  #filter(abs(value) < 500) %>%
  graph_from_data_frame(directed = FALSE)


graph_cors2
str(graph_cors2)
names(graph_cors2)
names(mtcars)
head(graph_cors, 30)
tail(graph_cors)
col <- c(rep("red", 5), rep("green", 6))

ggraph(graph_cors) +
  geom_edge_link(aes(edge_alpha = abs(r), 
                     edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), 
                              colors = c("firebrick2", "dodgerblue2")) +
  geom_node_point(color = col, size = 10) +
  geom_node_text(aes(label = name)) +
  theme_graph() +
  labs(title = "Correlations between car variables")

ggraph(graph_cors2) +
  geom_edge_link(aes(edge_alpha = abs(value/max(value)), 
                     edge_width = abs(value/max(value)), color = "blue")) +
  guides(edge_alpha = "none", edge_width = "none") +
  #scale_edge_colour_gradientn(limits = c(-1, 1), 
  #                            colors = c("firebrick2", "dodgerblue2")) +
  geom_node_point(color = col[1:11], size = 10) +
  geom_node_text(aes(label = name)) +
  theme_graph() +
  labs(title = "Hierarchical distances between car variables")


x1 <- rnorm(100)
y1 <- x1 + rnorm(100)
z1 <- x1 + rnorm(100)
x2 <- rnorm(100)
y2 <- x2 + rnorm(100)
z2 <- x2 + rnorm(100)
yo <- data.frame(x1, y1, z1, x2, y2, z2)
plot(yo)
mydata <- t(scale(yo))
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:5) wss[i] <- sum(kmeans(mydata, 
  	centers=i)$withinss)
plot(1:5, wss, type="b", xlab="Number of Clusters",
  ylab="Within groups sum of squares")

d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D2") 
plot(fit) # display dendogram
groups <- cutree(fit, k=2) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=2, border="red")

d
hist(d)
hist(cophenetic(fit))
d4a <- cophenetic(fit)
d4 <- myFun(d4a)
head(d4)
d42 <- 
d4<-d4[,c(1,2,3)]
hist(d4$value)
hist(log(d4$value,10))
d4$value <- log(d4$value,10)
graph_cors2 <- d4 %>%
  filter(abs(value) < 1.1) %>%
  graph_from_data_frame(directed = FALSE)

V(graph_cors2)$name

graph_cors3 <- t(mydata)  %>% correlate() %>% 
  stretch() %>%
  filter(abs(r) > 0.3) %>%
  graph_from_data_frame(directed = FALSE)


#d4 <- yo %>% 
#  correlate() %>% 
#  stretch() %>%
#  graph_from_data_frame(directed = FALSE)
#str(d4)

ggraph(graph_cors2) +
  geom_edge_link(aes(edge_alpha = 0.5, 
                     edge_width = 0.5, color = "blue")) +
  guides(edge_alpha = "none", edge_width = "none") +
  #scale_edge_colour_gradientn(limits = c(-1, 1), 
  #                            colors = c("firebrick2", "dodgerblue2")) +
  geom_node_point(color = "red", size = 10) +
  geom_node_text(aes(label = name)) +
  theme_graph() +
  labs(title = "Hierarchical distances between car variables")


ggraph(graph_cors3) +
  geom_edge_link(aes(edge_alpha = abs(r), 
                     edge_width = abs(r), color = "blue")) +
  guides(edge_alpha = "none", edge_width = "none") +
  #scale_edge_colour_gradientn(limits = c(-1, 1), 
  #                            colors = c("firebrick2", "dodgerblue2")) +
  geom_node_point(color = "red", size = 10) +
  geom_node_text(aes(label = name)) +
  theme_graph() +
  labs(title = "Correlations between car variables")

