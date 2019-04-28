
data = read.csv("womein-in-STEM.csv", header=TRUE, sep=",")  
# data from fivethirtyeight.com (https://github.com/fivethirtyeight/data)

attributes(data)

data$Median = NULL
data$Major_code = NULL
data$ShareWomen = NULL
data$Men = NULL
data$Rank = NULL

data$Major = as.numeric(factor(data$Major))
data$Major_category = as.numeric(factor(data$Major_category))

mat = data.matrix(data)
mat[c(1:5,70:75),]