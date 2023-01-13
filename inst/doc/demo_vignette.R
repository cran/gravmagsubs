## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      warning = FALSE,
                      fig.show = "hide",
                      fig.path = "figures/demo-")

## ----prereqs, eval=FALSE------------------------------------------------------
#  install.packages("fields")
#  install.packages("ggplot2")
#  install.packages("gridExtra")
#  install.packages("scales")
#  install.packages("scatterplot3d")

## ----import, results="hide"---------------------------------------------------
library(gravmagsubs)

## ----grvargs, eval=FALSE, results="hide"--------------------------------------
#  args(rectprismgrav)

## ----gravhelp, eval=FALSE, results="hide"-------------------------------------
#  help(rectprismgrav)

## ----gravquestion, eval=FALSE, results="hide"---------------------------------
#  ?rectprismgrav

## ----gravstation--------------------------------------------------------------
# location of the point where the gravity anomaly will be calculated
gravstation <- data.frame(x=0, y=0, z=0)

## ----prism1-------------------------------------------------------------------
# the rectangular prism is defined by its six edges
prism1 <- data.frame(xmin=-5, xmax=5,
                     ymin=-5, ymax=5,
                     zmin=-10, zmax=-5)

## ----dhro---------------------------------------------------------------------
# density contrast in g/cc
drho <- 0.3

## ----gravanom-----------------------------------------------------------------
gravanom <- rectprismgrav(gravstation$x, gravstation$y, gravstation$z,
                          prism1$xmin, prism1$xmax,
                          prism1$ymin, prism1$ymax,
                          prism1$zmin, prism1$zmax, drho)

## ----gravanom.echo------------------------------------------------------------
gravanom

## ----gravanom.redo------------------------------------------------------------
gravstation <- data.frame(x=-2, y=2, z=0)

gravanom <- rectprismgrav(gravstation$x, gravstation$y, gravstation$z,
                          prism1$xmin, prism1$xmax,
                          prism1$ymin, prism1$ymax,
                          prism1$zmin, prism1$zmax, drho)

## ----gravanom.redo.echo-------------------------------------------------------
gravanom

## ----grav.calc----------------------------------------------------------------
grav.calc <- data.frame(X = seq(-25, by=.1, length=500)) 

## ----stationsYZ---------------------------------------------------------------
grav.calc$Y <- 0
grav.calc$Z <- 0

## ----Geq----------------------------------------------------------------------
rho <- -0.67    # density contrast (g/cc)
gamma <- 6.674  # gravitational constant (x 1.e-11 m^3 / (kg s^2))
h1 <- 1.5       # infinite slab thickness (km)
h2 <- 1         # fault offset (km)
grav.calc$Geq <- 2*pi*rho*gamma*h1 + 2*rho*gamma*h2*(pi/2 + atan(grav.calc$X/2))

## ----Xzr----------------------------------------------------------------------
prism.width.h <- 1       # horizontal prism width (km)
prism.width.v <- 0.1     # vertical prism thickness (km)
X1 <- seq(-499.5, 499.5, by=prism.width.h)
Z1 <- seq(-3.95, -0.05, by=prism.width.v)
XZr <- expand.grid(xcenter=X1, Y=0, zcenter=Z1)

## ----plot.section067, message=FALSE, fig.dim = c(6,2)-------------------------
XZr$density <- 0
XZr$density[XZr$zcenter > -h1] <- rho
XZr$density[XZr$zcenter > -(h1 + h2) & XZr$xcenter > 0] <- rho

library(ggplot2)
library(scales)

# plot with custom color scale 
ggplot() +
  geom_raster(data = XZr, aes(x = xcenter, y = zcenter, fill = density)) +
  scale_fill_gradientn(colours=c("yellow", "darkred"),
                       values = rescale(c(-0.67, -0.33, 0)),
                       breaks = c(-0.67, 0),
                       labels = paste(c(-0.67, 0)),
                       limits = c(-0.67, 0),
                       name = "") +
  geom_point(data = grav.calc, aes(x=X, y=Z)) +
  labs(x = "Distance [km]",
       y = "Depth [km]",
       title = "Density contrast (g/cc))") +
  theme(panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size=16),
        aspect.ratio = 1/5)

## ----Ggms---------------------------------------------------------------------
grav.calc$Ggms <- as.vector(rectprismgrav(grav.calc$X, grav.calc$Y, grav.calc$Z,
                                          XZr$xcenter - prism.width.h/2,
                                          XZr$xcenter + prism.width.h/2,
                                          XZr$Y - 500, XZr$Y + 500,
                                          XZr$zcenter - prism.width.v/2,
                                          XZr$zcenter + prism.width.v/2, XZr$density) )

## ----Ggms.summary-------------------------------------------------------------
summary(grav.calc$Ggms)

## ----Geq.summary--------------------------------------------------------------
summary(grav.calc$Geq)

## ----sd.diff------------------------------------------------------------------
sd(grav.calc$Geq - grav.calc$Ggms)

## ----mean.diff----------------------------------------------------------------
mean(grav.calc$Geq - grav.calc$Ggms)

## ----plot.grav.calc, fig.dim=c(7,4)-------------------------------------------
ggplot() +
  geom_line(data = grav.calc, aes(x = X, y = Geq,
                                  colour = "Calculated from equation")) +
  geom_line(data = grav.calc, aes(x = X, y = Ggms,
                                  colour = "Calculated from rectprismgrav()")) +
  scale_colour_manual("", values=c("blue", "darkgreen")) +
  labs(title = "", x = "Distance [km]", y = "Gravity anomaly (mGal)") +
  theme(panel.background=element_rect(colour="black", fill="white"),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.2, 0.3),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold",hjust = 0.5, size=16 ),
        aspect.ratio = 1/4)

## ----density1-----------------------------------------------------------------
XZr$density1 <- 1

## ----sensmat------------------------------------------------------------------
XZr.sensmat <- rectprismgrav(grav.calc$X, grav.calc$Y, grav.calc$Z,
                             XZr$xcenter - prism.width.h/2,
                             XZr$xcenter + prism.width.h/2,
                             XZr$Y - 500, XZr$Y + 500,
                             XZr$zcenter - prism.width.v/2,
                             XZr$zcenter + prism.width.v/2, XZr$density1,
                             bycell=TRUE)

dim(XZr.sensmat)

## ----plot.sensmat, fig.dim=c(6,2)---------------------------------------------
# create a data frame with locations and log of the 10th gravity station
xzs1 <- data.frame(x=XZr[,1], z=XZr[,3], sensiv=log10(XZr.sensmat[10,]))

# create the plot
ggplot() +
  geom_raster(data = as.data.frame(xzs1) , aes(x = x, y = z, fill = sensiv)) +
  scale_fill_gradientn(colours=c("yellow","orange","red", "darkred"),
                       name="") +
  geom_point(data = grav.calc, aes(x=X, y=Z), colour="darkgray") +
  geom_point(data = grav.calc[10,], aes(x=X, y=Z), colour="black") +
  labs(x = "Distance [km]", 
       y = "Depth [km]",
       title = "log(Sensitivity(mGal per g/cc))") +
  theme(panel.background=element_rect(colour="black", fill="white"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold",hjust = 0.5, size=16 ),
        aspect.ratio=1/4)

## ----density2-----------------------------------------------------------------
XZr$density2 <- 0

## ----density3-----------------------------------------------------------------
XZr$density2[XZr$density == rho] <- -0.9
XZr$density3 <- 0
XZr$density3[XZr$density == rho] <- -0.4

## ----gravanom.models----------------------------------------------------------
gravanom.models <- XZr.sensmat %*% as.matrix(XZr[,c("density", "density2", "density3")])

## ----plot.gravanom.models, fig.dim=c(7,4)-------------------------------------
grav.calc$gmod.orig <- gravanom.models[,1]
grav.calc$gmod.low <- gravanom.models[,2]
grav.calc$gmod.high <- gravanom.models[,3]

ggplot() +
  geom_line(data = grav.calc, linewidth=2,
            aes(x = X, y = gmod.low, colour = "-0.9 g/cc")) +
  geom_line(data = grav.calc, linewidth=2,
            aes(x = X, y = gmod.orig, colour = "-0.67 g/cc")) +
  geom_line(data = grav.calc, linewidth=2,
            aes(x = X, y = gmod.high, colour = "-0.4 g/cc")) +
  scale_colour_manual(name="Density contrast", values=c("blue", "black", "red")) +
  labs(title = "", x = "Distance [km]", y = "Gravity anomaly (mGal)") +
  theme(panel.background=element_rect(colour="black", fill="white"),
        legend.key = element_rect(fill = "white"),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size=16 ),
        aspect.ratio = 1/4)

## ----fields, results="hide", message=FALSE------------------------------------
library(fields)
set.seed(123)

## ----rf.model, results="hide"-------------------------------------------------
# setting the anisotropy matrix
aniso.mat <- matrix(c(15/prism.width.h, 0, 0, 0.3/prism.width.v), nrow=2, byrow=TRUE)
# need to re-format the simulation grid using indices
rf.grid <- list(x = 1:length(X1), y = 1:length(Z1))
# setting up the random field model
rf.model <- stationary.image.cov(setup=TRUE, grid = rf.grid, V=aniso.mat )

## ----grav.sims, results="hide"------------------------------------------------
std.dev <- 0.14
grav.sim.1 <- sim.rf(rf.model) * std.dev
grav.sim.2 <- sim.rf(rf.model) * std.dev

## ----plot.random.sims, fig.dim=c(6,7)-----------------------------------------
library(gridExtra)
	
# Create the plot with the title: sim1 density contrast (g/cc)

xzg <- data.frame(x=XZr[,1], z=XZr[,3],
                  gravsim1=c(grav.sim.1),
                  gravsim2=c(grav.sim.2))

g1.plot<- ggplot() +
  geom_raster(data = as.data.frame(xzg), aes(x = x, y = z, fill = gravsim1)) +
  scale_fill_gradientn(colours=c("yellow", "orange", "red", "darkred"),
                       breaks = c(-0.5, 0, 0.5),
                       labels = paste(c(-0.5, 0 ,0.5)),
                       name = "") +
  geom_point(data = grav.calc, aes(x=X, y=Z)) +
  labs(x = "Distance [km]",
       y = "Depth [km]",
       title = "sim1 density contrast (g/cc)") +
  theme(panel.background=element_rect(colour="black", fill="white"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size=16 ),
        legend.position = "bottom",
        aspect.ratio = 1/5)
		
# Create the plot with the title: sim2 density contrast (g/cc)

g2.plot<- ggplot() +
  geom_raster(data = as.data.frame(xzg), aes(x = x, y = z, fill = gravsim2)) +
  scale_fill_gradientn(colours=c("yellow", "orange", "red", "darkred"),
                       breaks = c(-0.5, 0, 0.5),
                       labels = paste(c(-0.5, 0, 0.5)),
                       name = "") +
  geom_point(data = grav.calc, aes(x=X, y=Z)) +
  labs(x = "Distance [km]",
       y = "Depth [km]",
       title = "sim2 density contrast (g/cc)") +
  theme(panel.background = element_rect(colour="black", fill="white"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size=16 ),
        legend.position = "bottom",
        aspect.ratio = 1/5)
		
grid.arrange(g1.plot, g2.plot, nrow=2)


## ----grav.sims.anom-----------------------------------------------------------
grav.sims.mat <- matrix(c(grav.sim.1, grav.sim.2) , byrow=FALSE, ncol =2)
grav.sims.anom <- XZr.sensmat %*% grav.sims.mat
dim(grav.sims.anom)

## ----plot.grav.profile, fig.dim=c(6,6)----------------------------------------
anom.df <- data.frame(x = grav.calc$X,
                      anom1 = grav.sims.anom[,1],
                      anom2 = grav.sims.anom[,2])

g1.plot<- ggplot() +
  geom_line(data = as.data.frame(anom.df), aes(x = x, y = anom1)) +
  labs(x = "Distance [km]",
       y = "Anomaly (mGal)",
       title = "sim1 anomaly profile (mGal)") +
  theme(panel.background = element_rect(colour="black", fill="white"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size=16 ),
        legend.position = "bottom",
        aspect.ratio = 1/5)

g2.plot<- ggplot() +
  geom_line(data = as.data.frame(anom.df), aes(x = x, y = anom2)) +
  labs(x = "Distance [km]",
       y = "Anomaly (mGal)",
       title = "sim2 anomaly profile (mGal)") +
  theme(panel.background = element_rect(colour="black", fill="white"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size=16 ),
        legend.position = "bottom",
        aspect.ratio = 1/5)

grid.arrange(g1.plot, g2.plot, nrow=2)

## ----D3.stations--------------------------------------------------------------
xstations <- seq(-5, 5, by=0.1)
ystations <- seq(-5, 5, by=0.1)
D3.stations <- expand.grid(xstations, ystations)
D3.stations$Zstation <- 0
names(D3.stations) <- c("Xstation", "Ystation", "Zstation")

## ----plot.scatter3d, fig.dim=c(6,5)-------------------------------------------
width <- 0.1
half.width <- width / 2

xsource <- seq(-2 + half.width, 2 - half.width, by=width)
ysource <- seq(-2 + half.width, 2 - half.width, by=width)
zsource <- seq(-3 + half.width, -1 - half.width, by=width)

D3.source <- expand.grid(xsource, ysource, zsource)
names(D3.source) <- c("xcenter", "ycenter", "zcenter")

D3.source$density <- 0.1 * D3.source$ycenter

# define the color ramp by splitting the density model into 10 equally-spaced intervals
yorPal <- colorRampPalette(c('khaki', 'orange', 'firebrick'))
D3.cols = cut(D3.source$density, breaks=10)

img3d <- scatterplot3d::scatterplot3d(D3.source$xcenter,
                                      D3.source$ycenter,
                                      D3.source$zcenter,
                                      pch=16, cex.symbols=.5,
                                      color=yorPal(10)[as.numeric(D3.cols)],
                                      xlab="Distance [km]", ylab="[km]", zlab="[km]")

leg.minmax <- c(paste("+", max(D3.source$density), " g/cc", sep=""),
                " 0.0 g/cc",
                paste(min(D3.source$density), "g/cc"))
legend(img3d$xyz.convert(3.7, -2, -2.7), legend = leg.minmax,
       col=c("firebrick", "orange", "khaki"), pch=16, title="Density contrast",
       inset = -0.25, xpd = TRUE, bty="n")

## ----plot.grav3d, fig.dim=c(5,5)----------------------------------------------
D3.gravanom <- rectprismgrav(D3.stations$Xstation,
                             D3.stations$Ystation,
                             D3.stations$Zstation,
                             D3.source$xcenter - half.width,
                             D3.source$xcenter + half.width,
                             D3.source$ycenter - half.width,
                             D3.source$ycenter + half.width,
                             D3.source$zcenter - half.width,
                             D3.source$zcenter + half.width,
                             D3.source$density,
                             bycell=FALSE)

# create data frame
d3.st <- data.frame(x=D3.stations$Xstation, y=D3.stations$Ystation,
                    val1 = D3.gravanom)

# plot with custom color scale 
ggplot() +
  geom_raster(data = d3.st, aes(x = x, y = y, fill = val1)) +
  scale_fill_gradientn(
    colours = c("darkblue", "blue", "skyblue", "white", "tomato", "red", "darkred"),
    values = rescale(c(-1.5, 0, 1.5)),
    breaks = seq(-1.5, 1.5, 0.5),
    labels = paste(seq(-1.5, 1.5, 0.5)),
    limits = c(-1.6, 1.6),
    name = "") +
  labs(x = "Easting [km]", 
       y = "Northing [km]",
       title = "Gravity anomaly (mGal)") +
  theme(panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size=16),
        aspect.ratio = 1)

## ----magargs, eval=FALSE, results="hide"--------------------------------------
#  args(rectprismmag)

## ----maghelp, eval=FALSE, results="hide"--------------------------------------
#  help(rectprismmag)

## ----magquestion, eval=FALSE, results="hide"----------------------------------
#  ?rectprismmag

## ----D3.suscnorm--------------------------------------------------------------
D3.source$suscnorm <- rnorm(n=length(D3.source[,1]), mean=0.015, sd=0.001)

## ----plot.mag3d.ind, fig.dim=c(5,5)-------------------------------------------
D3.source$nrmstr <- 0
D3.source$nrmdecl <- 0
D3.source$nrmincl <- 0
D3.source$fieldtotal <- 48800
D3.source$fielddecl <- 12
D3.source$fieldincl <- 60

D3.maganomind <- rectprismmag(D3.stations$Xstation,
                              D3.stations$Ystation,
                              D3.stations$Zstation,
                              D3.source$xcenter - half.width,
                              D3.source$xcenter + half.width,
                              D3.source$ycenter - half.width,
                              D3.source$ycenter + half.width,
                              D3.source$zcenter - half.width,
                              D3.source$zcenter + half.width,
                              suscvolsi = D3.source$suscnorm,
                              nrmstr = D3.source$nrmstr,
                              nrmdecl = D3.source$nrmdecl,
                              nrmincl = D3.source$nrmincl,
                              fieldtotal = D3.source$fieldtotal,
                              fielddecl = D3.source$fielddecl,
                              fieldincl = D3.source$fieldincl,
                              bycell=FALSE)

# create data frame
d3.st <- D3.stations
d3.st$val <- D3.maganomind

# plot with custom color scale
ggplot() +
  geom_raster(data = as.data.frame(d3.st),
              aes(x = Xstation, y = Ystation, fill = val)) +
  scale_fill_gradientn(
    colours = c("darkblue", "blue", "lightblue", "white", "pink", "red", "darkred"),
    values = rescale(c(min(d3.st$val), 0, max(d3.st$val))),
    breaks = c(-50, 0, 50, 100),
    labels = paste(c("-50", 0, "50", "100")),
    limits = c(-50, max(d3.st$val)),
    name = "") +
  labs(x = "Distance [km]", 
       y = "Distance [km]",
       title = "Induced magnetic anomaly (nT)") +
  theme(panel.background = element_rect(colour="black", fill="white"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size=16),
        aspect.ratio = 1)


## ----plot.mag3d.nrm.gonzalez, fig.dim=c(5,5)----------------------------------

xstations <- seq(-5, 10, by=0.25)
ystations <- seq(-6, 10, by=0.25)
stations <- expand.grid(xstations, ystations)
stations$Zstation <- 0
names(stations) <- c("Xstation", "Ystation", "Zstation")

maganomnrm <- rectprismmag(stations$Xstation,
                           stations$Ystation,
                           stations$Zstation,
                           xmin = 1,
                           xmax = 3,
                           ymin = 1,
                           ymax = 3.5,
                           zdeep = -2.5,
                           zshallow = -0.55,
                           suscvolsi = 0.1,
                           nrmstr = 5,
                           nrmdecl = -15,
                           nrmincl = -30,
                           fieldtotal = 23639,
                           fielddecl = -30,
                           fieldincl = -45,
                           bycell=FALSE)

# create data frame
d3.st <- data.frame(x=stations[,1], y=stations[,2], val1 = maganomnrm)

# plot with custom color scale 
ggplot() +
  geom_raster(data = d3.st, aes(x = x, y = y, fill = val1)) +
  scale_fill_gradientn(
    colours = c("darkblue", "blue", "skyblue", "white", "tomato", "red", "darkred"),
    values = rescale(c(-1000, 0,1500)),
    breaks = c(-1000, -500, 0, 500, 1000, 1500),
    labels = paste(c(-1000, -500, 0, 500, 1000, 1500)),
    limits = c(-1000, 1500),
    name = "") +
  labs(x = "Easting [km]", 
       y = "Northing [km]",
       title = "Total magnetic anomaly (nT)") +
  theme(panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size=16),
        aspect.ratio = 1)


## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

