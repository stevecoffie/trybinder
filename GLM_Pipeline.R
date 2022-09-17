library(GLMsData)
library(faraway)
library(statmod)
library(ggplot2)
library(dplyr)
library(GGally)
library(corrplot)
library(Hmisc)
library(reshape)

# importing data
data(pipeline)
pip <- pipeline
attach(pip)
head(pip, 2)

# missing variables
sum(is.na(pip))

# dimension
dim(pip)

# variable types
str(pip)

# summary
describe(pip)

# scatter plot
ggpairs(pip[,1:2], lower=list(continuous="smooth"))+ theme_bw() + 
  labs(title="Cancer Mean") + 
  theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=12))

# Basic histogram
ggplot(pip, aes(x=Lab)) + geom_histogram(binwidth=15)
ggplot(pip, aes(x=Field)) + geom_histogram(binwidth=15)


ggplot(data = pip, aes(x=Batch, y=Lab)) + 
  geom_boxplot(aes(fill=Batch)) + 
  xlab("Variables") + ylab("")+ guides(fill=guide_legend(title="Batch"))

# scatter plot
plot( log(Lab) ~ Field, las=1, pch=19,
      xlab="Field", ylab="Lab")


# Now compute means and variances of each origin/age group:
pip$Field1 <- cut(pip$Field, breaks = 5)
mn <- with( pip, tapply(Lab, Field1, mean ) )
mn
vr <- with( pip, tapply(Lab, Field1, var ) )
vr
# Plot
plot( vr ~ mn, las=1, pch=19,
      xlab="group means", ylab="roup variance")
# log plot
plot( log(vr) ~ log(mn), las=1, pch=19,
      xlab="log(group means)", ylab="log(group variance)")

plot( sqrt(vr) ~ mn, las=1, pch=19,
      xlab="log(group means)", ylab="log(group variance)")

data.frame(mn, vr, ratio=vr/mn)

mvline <- lm( c(vr)  ~  c(mn) ) 
slope <- round( coef( mvline )[2], 2)
abline( mvline, lwd=2)
slope


# model
pip.mod <- glm(Lab ~ Field, family = Gamma(link = "log"), data = pip)
summary(pip.mod)

# deviance residual distribution
#deviance residual is approximate normal,  phi less than 1/3
# dispersion parameter
summary(pip.mod)$dispersion > 1/3

# standardized deviance residuals against the fitted values
scatter.smooth( rstandard(pip.mod) ~ fitted(pip.mod), las=1, 
                ylab="Standardized deviance residual", 
                xlab="Fitted values" )
# standardized deviance residuals against  fitted values transformed
# to the constant-information scale
scatter.smooth( rstandard(pip.mod) ~ log(fitted(pip.mod)), las=1, 
                ylab="Standardized deviance residual", 
                xlab="log(Fitted values)" )

# working residuals
z <- resid(pip.mod, type="working") + pip.mod$linear.predictor
scatter.smooth( z ~ pip.mod$linear.predictor, las=1, 
                xlab="Working responses, z", ylab="Linear predictor")
abline(0, 1, col="blue") # Adds line of equality

# partial residual
termplot(pip.mod, partial.resid=TRUE, las=1)

# Q-Q Plot
qr.pip <- qresid( pip.mod )
qqnorm( qr.pip, las=1 )
qqline( qr.pip)

# cook distance
pip.cd <- cooks.distance( pip.mod)
plot( pip.cd, type="h", ylab="Cook's distance", las=1)














lm.fit <- lm(Lab ~ Field, data = pipeline)
summary(lm.fit)
plot(pipeline$Field, pipeline$Lab, col=as.factor(pipeline$Batch),main = ("Lab ~ Field"),pch='*', cex=1.5)
legend("topright",title="Batch",legend= unique(pipeline$Batch), fill=1:length(pipeline$Batch) )
abline(coef(lm.fit),lty=5)


plot(residuals(lm.fit) ~ Field, pipeline, main ="Residuals versus log(time) for simple linear model")

i <- order(pipeline$Field) 
npipe <- pipeline[i,] 
ff <- gl(12,9)[-108] 
meanfield <- unlist(lapply(split(npipe$Field,ff),mean))
varlab <- unlist(lapply(split(npipe$Lab,ff),var))

plot(log10(varlab) ,log10(meanfield))
df<-data.frame(cbind(as.numeric(meanfield),as.numeric(varlab)))
df <- df[-c(12),]
lm.var.model <- lm(log(varlab) ~ log(meanfield) , data= df)
summary(lm.var.model)

