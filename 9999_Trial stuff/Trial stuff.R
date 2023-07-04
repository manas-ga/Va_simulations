library(ggplot2)
d = read.table("Mutations.txt", header=T, sep = " ")
d$Frequency = d$Number/20000
d$Diversity = 2*d$Frequency*(1-d$Frequency)

p = ggplot(d, aes(x = Diversity, y = s*s))
p + theme_bw() + geom_point() + geom_smooth(method = "lm")

hist(d$Frequency, breaks = 20)
hist(d$Diversity, breaks = 20)
hist(d$s, breaks = 20)

fit_1 = lm(s*s ~ Diversity, d)
summary(fit_1)
