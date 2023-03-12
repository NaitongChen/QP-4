library(ggplot2)

df = read.csv("../data/result.csv", sep = ",")

colnames(df)[1] = "Region"

ggplot(df, aes(x=Region, y=PointEstimate)) + 
geom_pointrange(aes(ymin=lb, ymax=up))