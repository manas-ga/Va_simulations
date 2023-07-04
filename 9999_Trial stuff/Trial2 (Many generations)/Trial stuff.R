#library(ggplot2)
d = read.table("mutations_all_gen.txt", header=T, sep = " ")
d$Frequency = d$Number/20000
d$Diversity = 2*d$Frequency*(1-d$Frequency)
d$V_a = (d$Diversity)*(d$s)*(d$s)


# Create an empty data frame to store the outputs of linear models between log(s*s) and log(Diversity) for each generation


alpha_div = data.frame()




# List various generations in the data set

gen_list = levels(factor(d$Tick_Number))

#loop over generations

for (generation in gen_list){

fit_1 = lm(log(s*s) ~ log(Diversity), d[d$Tick_Number==generation,])
#print(generation)
#summary(fit_1)

alpha_div = rbind(alpha_div,c(generation, coef(fit_1)[1], coef(fit_1)[2], length(d[d$Tick_Number==generation,]$Tick_Number), mean(d[d$Tick_Number==generation,]$Diversity), sum(d[d$Tick_Number==generation,]$V_a)))


}


colnames(alpha_div) = c("Generation", "Intercept", "Slope", "Segregating_mutations", "Average_diversity", "V_a")


pdf("r_4e-7.pdf")
plot(alpha_div$Generation, alpha_div$Slope)
plot(alpha_div$Generation, alpha_div$Intercept)
plot(alpha_div$Generation, alpha_div$Segregating_mutations)
plot(alpha_div$Generation, alpha_div$Average_diversity)
plot(alpha_div$Generation, alpha_div$V_a)
plot(alpha_div$Segregating_mutations, alpha_div$Average_diversity)
plot(alpha_div$Slope, alpha_div$Intercept)
dev.off()
#############################################

#p = ggplot(d, aes(x = Diversity, y = s*s))
#p + theme_bw() + geom_point() + geom_smooth(method = "lm")

#p = ggplot(d, aes(x = log(Diversity), y = log(s*s))
#p + theme_bw() + geom_point() + geom_smooth(method = "lm")

#hist(d$Frequency, breaks = 20)
#hist(d$Diversity, breaks = 20)
#hist(d$s, breaks = 20)
