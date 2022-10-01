# Example analyses using qte_by_outcome() and plot_qte_by_outcome() functions.
# Harri Hemila, Matti Pirinen
# Oct 1, 2022

# Read in qte_by_outcome functions by source(filename) command
# where filename is path to the file 'qte_by_outcome_functions.R'
# in your computer.
# Alternatively, you can read the file directly from GitHub by using
filename = "https://raw.githubusercontent.com/mjpirinen/qte-by-outcome/main/qte_by_outcome_functions.R" 
source(filename)

# Data files can be found from GitHub 
# https://github.com/mjpirinen/qte-by-outcome
# Below direct web links are provided.



# Analysis of Mossad data on effect of zinc gluconate lozenges 
# on duration of common cold.
# Data taken from pages 2-4 of the additional file 2 of
# https://doi.org/10.1186/s12874-017-0356-y

filename.1 = "https://raw.githubusercontent.com/mjpirinen/qte-by-outcome/main/Mossad.csv"
Mossad <- read.csv(filename.1)
str(Mossad) # Check data looks OK.
table(Mossad$Duration, Mossad$Zinc)

Treatment <- subset(Mossad, Zinc == 1)$Duration
table(Treatment)
Control <- subset(Mossad, Zinc == 0)$Duration
table(Control)

set.seed(1) #Set random seed for reproducible results
res.qq = qte_by_outcome(Treatment = Treatment, Control = Control) #Estimate treatment effect

#Check output of qte_by_outcome() function
res.qq

par(mfrow = c(1,2))

#Plot results on outcome scale ('plot.rte = F')
plot_qte_by_outcome(res.qq, plot.rte = F, 
                    xlim = c(0,20), xaxs = "i", col = "red",
                    ylab = "Treatment effect (Days)", 
                    xlab="Duration in the control group (Days)")
abline(-4.0, 0, col = "blue", lt = 3,  lw = 2) #Add known mean effect line
text(18, 0, "A", pos = 1, cex = 3) #Label this one as panel "A"

#Plot results on relative scale ('plot.rte = T')
plot_qte_by_outcome(res.qq, plot.rte = T, 
                    xlim = c(0,20), xaxs = "i", col = "red",
                    xlab="Duration in the control group (Days)")
abline(-43, 0, col = "blue", lt = 3, lw = 2) #Add known mean effect line
text(18, 0, "B", pos = 1, cex = 3) #Label this one as panel "B"


#
##
###
##
#


# Analysis of data of the three randomized trials 
# on zinc acetate lozenges on duration of common cold.
# The data are from page 8 of supplementary file 2 of
# https://doi.org/10.1093/ofid/ofx059

filename.2 = "https://raw.githubusercontent.com/mjpirinen/qte-by-outcome/main/ZnAcet.csv"
ZnAcet <- read.csv(filename.2)
str(ZnAcet)
table(ZnAcet$Duration, ZnAcet$Zinc)

Treatment <- subset(ZnAcet, Zinc == 1)$Duration
Control <- subset(ZnAcet, Zinc == 0)$Duration

set.seed(2)
res.qq = qte_by_outcome(Treatment = Treatment, Control = Control)


par(mfrow = c(1,2))
plot_qte_by_outcome(res.qq, plot.rte = F,
                    xaxs = "i", xlim = c(0,20), col = "red", 
                    ylab = "Treatment effect (Days)", 
                    xlab="Duration in the control group (Days)")
abline(-2.7, 0, col = "blue", lt = 3,  lw = 2)
text(18, 0, "A", pos = 1, cex = 3)

plot_qte_by_outcome(res.qq, plot.rte = T,
                    xaxs = "i", xlim = c(0,20), col = "red", 
                    xlab="Duration in the control group (Days)")
abline(-36, 0, col = "blue", lt = 3, lw = 2)
text(18, 0, "B", pos = 1, cex = 3)


#
##
###
##
#


# Analysis of data of the two nasal iota-carrageenan trials 
# on duration of common cold.
# The data are from page 3 of the supplementary file of
# https://doi.org/10.1002/prp2.810
# sensored data on day 20 are imputed with day 20 
# (6 in carrageenan, 21 in placebo).


filename.3 = "https://raw.githubusercontent.com/mjpirinen/qte-by-outcome/main/Carrageenan.csv"
Carrageenan <- read.csv(filename.3)
str(Carrageenan)
table(Carrageenan$Duration,Carrageenan$Carr)

Treatment <- subset(Carrageenan, Carr == 1)$Duration
Control <- subset(Carrageenan, Carr == 0)$Duration

set.seed(3)
res.qq = qte_by_outcome(Treatment = Treatment, Control = Control, 
                        at = c(3:19)) #at specifies the estimation points

par(mfrow = c(1,2))
plot_qte_by_outcome(res.qq, plot.rte = F, 
                    xlim = c(0,20), ylim = c(-9,1),
                    yaxp = c(-8, 0, 4),
                    col = "red", xaxs = "i", 
                    ylab = "Treatment effect (Days)", 
                    xlab="Duration in the control group (Days)")
abline(0, 0, col = "black", lt = 2, lw = 1)
abline(-1.9, 0, col = "blue", lt = 3, lw = 2)
text(18, 0, "A", pos = 1, cex = 3)

plot_qte_by_outcome(res.qq, plot.rte = T, 
                    xlim = c(0,20),  ylim = c(-60,25),
                    yaxp = c(-60, 20, 4), 
                    col = "red", xaxs = "i",
                    xlab="Duration in the control group (Days)")
abline(0, 0, col = "black", lt = 2, lw = 1)
text(18, 18, "B", pos = 1, cex = 3)



# Comparison of the Carrageenan data with imputation at 20 days 
# (continuous) against imputation at 30 days (dash line).
# Only difference to above carrageenan data is that 
# the censored observations on day 20 are now imputed to be 30 days


filename.4 = "https://raw.githubusercontent.com/mjpirinen/qte-by-outcome/main/Carrageenan30.csv"
Carrageenan30 <- read.csv(filename.4)
table(Carrageenan30$Duration, Carrageenan30$Carr)  
str(Carrageenan30)

Treatment30 <- subset(Carrageenan30, Carr == 1)$Duration
Control30 <- subset(Carrageenan30, Carr == 0)$Duration

set.seed(4)
#Use more bootstrap samples to ensure max accuracy (20000 while default = 2000)
res.qq = qte_by_outcome(Treatment = Treatment, Control = Control, 
                        B = 20000, at = c(3:19)) 
res30.qq = qte_by_outcome(Treatment = Treatment30, Control = Control30, 
                          B = 20000, at = c(3:19))
res30.qq

par(mfrow = c(1,1))
plot (res.qq$at, res.qq$qte, type ="l", xlim = c(0,20), 
      yaxp = c(-10, 0, 5), xaxs = "i",
      ylab = "Mean treatment effect (Days)", 
      xlab="Duration in the control group (Days)")
lines(res30.qq$at, res30.qq$qte, lty = 2)

# Conclusion: 
# No difference observed due to imputation in the range considered here.

