suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))

# Post-process the Q* file.
load('survival/data_dump/Q_star_ami.object')
# Q.star.list epsilon.mat tte.mat

head((table(round(sort(tte.mat[,1])))),10)
 
head((table(round(sort(tte.mat[,2])))),30)
head((table(round(sort(tte.mat[,3])))),30)
head((table(round(sort(tte.mat[,4])))),30)
head((table(round(sort(tte.mat[,5])))),30)

melted.tte <- melt(data.frame(tte.mat, check.names=FALSE), 
                   id.vars = NULL,
                   variable.name = "hospital",
                   value.name="tte")

table()

# All on top.
ggplot(data=melted.tte, aes(x=tte, colour=hospital)) + 
  geom_density()

ggplot(data=melted.tte, aes(x=tte, colour=hospital)) + 
  geom_density() + xlim(c(0,60))


?geom_density


ggplot(data=melted.tte, aes(x=tte)) + 
  facet_grid(hospital ~ .) +
  geom_density() 
  
ggplot(data=melted.tte, aes(x=tte,colour=hospital)) + 
  facet_grid(hospital ~ .) +
  xlim(c(0,60)) +
  geom_density() 

ggplot(data=melted.tte, aes(x=tte,colour=hospital)) + 
  facet_grid(hospital ~ .) +
  xlim(c(0,120)) +
  geom_density() 

ggplot(data=melted.tte, aes(x=tte,colour=hospital)) + 
  facet_grid(hospital ~ .) +
  geom_density() 


hist(tte.mat[,1])

round(tte.mat[,1][tte.mat[,1]<30])

