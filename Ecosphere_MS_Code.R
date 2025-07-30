# Code to reproduce analyses in Ecosphere manuscript:
# "Crossing boundaries: Introduced trout alter the bird community in a naturally fishless headwaters ecosystem"
# Authors: Mary K. Clapp, Erik W. Meyer, Gail L. Patricelli
# All code written by MKC


# LIBRARIES AND DATA ------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(lme4)
library(ggeffects)
library(ggsignif)
library(betapart)
library(glmmTMB)
library(chron)
library(vegan)
library(sjPlot)

# this function tests for overdispersion (the variance in the residuals are greater than the mean). if the ratio > 1, overdispersion 
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# data

mayfly <- read_csv("data/insects/mayflycounts.csv") 

aou <- read_csv("data/NACC_list_species.csv") 

ibp <- read_csv("data/IBP-AOS-LIST23.csv")

aoufull <- full_join(aou, ibp, by=c("common_name"="COMMONNAME")) %>% dplyr::select(common_name, order, family, genus, species, SPEC)

d3 <- read_csv("data/pointcount_data.csv")
lakes <- read_csv("data/SEKI_Lake_Metadata_2020.csv")
spec_names <- read_csv("data/species_names.csv") 

# format raw point count data
d3$Time <- times(d3$Time, format="hms")

dat <- d3 %>% 
  group_by(basin, fish, location, year, date, visit, point, wind, Time, jday) %>% 
  summarise(nspecies=n_distinct(species), count=n()) %>% 
  left_join(lakes, by = c("location"="LAKENAME")) %>%
  dplyr::select(basin, fish, location, year, jday, Time, count, wind, nspecies, ELEVATION) 
colnames(dat)

# scale independent variables
dat$fish <- factor(dat$fish, levels = c("fishless", "fish"))
dat$year <- factor(dat$year, levels = c("2020", "2015", "2014"))
dat$loc.pt <- paste(dat$location, dat$point, sep=".")

dat2 <- dat # dat2 will now be the dataset containing CENTER2R
dat <- dat2 %>% filter(location!="CENTER2R") # dat does not contain post-fish-removal bird surveys

# set themes and settings

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
fishcolors <- c(cbPalette[6], cbPalette[2]) # corresponds to fish as factor with levels (fishless, fish)

# ABUNDANCE GLMM ----------------------------------------------------------

summary(abun3 <- glmer(count ~ fish + scale(ELEVATION) +(1|year) + (1|basin/location/loc.pt), family = poisson, data = dat, glmerControl(calc.derivs = F, optCtrl=list(maxfun=10^4))) ) 

summary(abun0 <- glmer(count ~ scale(ELEVATION) +(1|year) + (1|basin/location/loc.pt), family = poisson, data = dat, glmerControl(calc.derivs = F, optCtrl=list(maxfun=10^4))) ) 

anova(abun0, abun3) 

# demonstrate that a poisson distribution is appropriate (corrects overdispersion)
abun_lm <- lmer(count ~ fish + scale(ELEVATION) +(1|year) + (1|basin/location/loc.pt), data = dat) 
overdisp_fun(abun_lm) # ratio 

# Examine Residuals
residuals <- residuals(abun3)
fitted_values <- fitted(abun3)

# 1. Plot residuals vs fitted values
plot(residuals ~ fitted_values)

# 2. Plot residuals vs each covariate
plot(dat$ELEVATION ~ residuals)
plot(dat$fish ~ residuals)

# create model table (Table S2)
tab_model(abun3, 
          transform = NULL,
          show.se = TRUE,
          dv.labels = c("model of count"),
          pred.labels = c("Intercept", "Fish", "Elevation"),
          string.se = "SE",
          string.ci = "95% CI",
          string.est = "Estimate",
          p.style = "numeric_stars",
          show.icc = F) 

# retrieve marginal means +/- 95% CI for fish effect
abun.fx <- ggpredict(abun3,terms = c("fish"), ci.level = 0.95)

dat %>% group_by(fish) %>% summarise(median=median(count))

abun.plot <- ggplot() +
  geom_jitter(width = 0.3, height = 0.2, data=dat, aes(x=fish, y=count, color=fish), alpha=0.3) +
  geom_violin(data=dat, aes(x=fish, y=count, color=fish), alpha=0.3) +
  geom_point(data=abun.fx, aes(x=x, y=predicted)) +
  geom_linerange(data=abun.fx, aes(x=x, ymin=conf.low, ymax=conf.high)) +
   geom_text(data = abun.fx, aes(x = x, y = predicted, label = round(predicted,2)),color="black", size = 4, vjust = 5) +
 ggpubr::geom_bracket(
 tip.length = 0.02, # the downard "tips" of the bracket
 vjust = 0, # moves text label (in this case, the p-value)
 xmin = 1, #starting point for the bracket
 xmax = 2, # ending point for the bracket
 y.position = 17.4, # vertical location of the bracket
 label.size = 4,# size of bracket text
 label = paste0("p = 0.037") ) +
  scale_color_manual(values=fishcolors) +
#  geom_text(aes(x=1.5, y=10, label = paste("p = ", format(p_value, digits=2)) )) +
  theme_pubr() +
  labs(x=NULL, y="number of detected individuals", colour = "Lake Type") +
  theme(text = element_text(size=12, family = "Arial"))
abun.plot

# species richness GLMM ---------------------------------------------------

summary(rich3 <- glmer(nspecies ~ fish + scale(ELEVATION) + (1|year) + (1|basin/location/loc.pt), family = poisson, data = dat, glmerControl(calc.derivs = F, optCtrl=list(maxfun=10^4))) ) 

summary(rich0 <- glmer(nspecies ~ scale(ELEVATION) + (1|year) + (1|basin/location/loc.pt), family = poisson, data = dat, glmerControl(calc.derivs = F, optCtrl=list(maxfun=10^4))) )  

anova(rich0,rich3)

tab_model(rich3, 
          transform = NULL,
          show.se = TRUE,
          dv.labels = c("model of count"),
          pred.labels = c("Intercept", "Fish", "Elevation"),
          string.se = "SE",
          string.ci = "95% CI",
          string.est = "Estimate",
          p.style = "numeric_stars",
          show.icc = F) #,
     #     file = "./tables/richMod_2025_4_Jun.doc") 

# retrieve marginal means +/- 95% CI for fish effect
rich.fx <- ggpredict(rich3 ,terms = c("fish"), ci.level = 0.95)


rich.plot <- ggplot() +
  geom_jitter(width = 0.3, height = 0.2, data=dat, aes(x=fish, y=nspecies, color=fish), alpha=0.3) +
  geom_violin(data=dat, aes(x=fish, y=nspecies, color=fish), alpha=0.3) +
  geom_point(data=rich.fx, aes(x=x, y=predicted)) +
  geom_linerange(data=rich.fx, aes(x=x, ymin=conf.low, ymax=conf.high)) +
   geom_text(data = rich.fx, aes(x = x, y = predicted, label = round(predicted,2)),color="black", size = 4, vjust = 5) +
  ggpubr::geom_bracket(
    tip.length = 0.02, # the downard "tips" of the bracket
    vjust = 0, # moves your text label (in this case, the p-value)
    xmin = 1, #starting point for the bracket
    xmax = 2, # ending point for the bracket
    y.position = 11, # vertical location of the bracket
    label.size = 4, # size of your bracket text
    label = paste0("p = 0.069") ) +
  scale_color_manual(values=fishcolors) +
  #  geom_text(aes(x=1.5, y=10, label = paste("p = ", format(p_value, digits=2)) )) +
  theme_pubr() +
  labs(x=NULL, y="number of detected species", colour = "Lake Type") +
  theme(text = element_text(size=12, family = "Arial"))
rich.plot

abunrich_tabs <- tab_model(abun3, rich3,
                           transform = NULL,
                           show.se = TRUE,
                           dv.labels = c("Abundance Model", "Richness Model"),
                           pred.labels = c("Intercept", "Fish", "Elevation"),
                           string.se = "SE",
                           string.ci = "95% CI",
                           string.est = "Estimate",
                           p.style = "numeric",
                           show.icc = F) #,
                       #    file = "./tables/abunrichMods_4Jun2025.doc") 
bothplots <- ggarrange(abun.plot, rich.plot,
                       legend = "none",
                       common.legend = T,
                       nrow=1,
                       labels = c("(a)","(b)"))

ggsave(plot=bothplots, filename = paste0("images/", "Figure2_", Sys.Date(), ".tiff"),
        width = 18, height = 12, units = "cm", dpi = 600)


# Beta Diversity ----------------------------------------------------------
# 2020 data only, removes trout-removal lake, summarizes by LAKE-level counts

# exact results may vary by run because set.seed() was not used.

nrow(d3 %>% filter(year==2020) %>% filter(location != "CENTER2R") %>% group_by(location, point) %>% summarise(nspec=n_distinct(species)))
# 75 total surveys

relabun <- d3 %>% filter(year==2020) %>% filter(location != "CENTER2R") %>% mutate(locpt = paste(location, point)) %>% group_by(species) %>% summarise(nsurv = n_distinct(locpt), propsurv = nsurv/75) %>% arrange(propsurv)


only20 <- d3 %>% filter(year==2020) %>% filter(location != "CENTER2R", species !="PIKA") %>%
  group_by(basin, location, fish, species) %>%
  summarise(lakecount = n()) %>%
  left_join(lakes, join_by(location==LAKENAME)) %>%
  mutate(relcount = lakecount/NPOINTS) %>%
 dplyr::select(-lakecount) %>%
  spread(key = species, value = relcount, fill = 0) 


mat <- only20 %>% as.data.frame() %>% dplyr::select(AMDI:YRWA)

n1 = metaMDS(mat, distance = "bray", trymax = 1000, k=3) 
stressplot(n1)
n1$stress
MDS1 = n1$points[,1]
MDS2 = n1$points[,2]
MDS3 = n1$points[,3]

# for plotting in ggplot
scores = data.frame(MDS1 = MDS1, MDS2 = MDS2, Fish = only20$fish, Basin = only20$basin, Lake = only20$location, Elev = only20$ELEVATION)

# species scores
species.scores <- as.data.frame(scores(n1, "species"))  #Uses the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)
species.scores <- left_join(species.scores, relabun)

# create spider plot
cent <- aggregate(cbind(MDS1, MDS2) ~ Fish, data = scores, FUN = mean)
segs <- merge(scores, setNames(cent, c('Fish','oNMDS1','oNMDS2')),
              by = 'Fish', sort = FALSE)

supercent <- cbind(mean(cent$MDS1), mean(cent$MDS2))

NMDSplot  <- ggplot(scores, aes(x = MDS1, y = MDS2)) +
  geom_segment(data = segs, aes(xend = oNMDS1, yend = oNMDS2, colour = Fish)) + 
  geom_point(data = cent, aes(colour = Fish, shape=Fish), size = 4) + 
  geom_point(aes(colour = Fish, shape = Fish)) +      
  geom_text(data=species.scores[species.scores$nsurv>2,],aes(x=NMDS1,y=NMDS2,label=species, alpha=propsurv), show.legend = F) +  # add the species labels
  scale_colour_manual(values=c("fish" = cbPalette[2], "fishless" = cbPalette[6])) +
  scale_shape(guide="none") +
   scale_x_continuous(breaks=seq(-1.25, 1.75, 0.5)) +
   scale_y_continuous(breaks=seq(-1.25, 1.75, 0.5)) +
  coord_equal() +
  theme_classic() +
  labs(colour="Lake Type", x="NMDS1", y="NMDS2", alpha=NULL) +
  theme_pubr() +
  theme(text = element_text(size=12, family = "Arial"), legend.position = "right",
        axis.text = element_text(size=10))

NMDSplot

#ggsave(NMDSplot, filename="images/forppt/NMDS_4June2025.png", width=6.5, height=5, dpi=600)
# Permanova ---------------------------------------------------------------

only20$fish <- factor(only20$fish, levels=c("fishless", "fish"))

bpart <- beta.pair.abund(mat, index.family = "bray")

adon.turn <- adonis2(bpart$beta.bray.bal ~ fish + scale(ELEVATION), data=only20)
adon.nest <- adonis2(bpart$beta.bray.gra ~ fish + scale(ELEVATION),  data=only20)
adon.bray <- adonis2(bpart$beta.bray ~ fish + scale(ELEVATION), data=only20)

adon.bray
adon.nest
adon.turn

scores$Fish <- factor(scores$Fish, levels=c("fishless", "fish"))
summary(lm(MDS1 ~  Elev, data=scores))
resids <- residuals(lm(MDS1 ~ Elev, data=scores))
summary(lm(resids ~ scores$Fish))

summary(nmds1.mod <- lm(MDS1 ~ Fish + scale(Elev), data=scores))
summary(nmds2.mod <- lm(MDS2 ~ Fish + scale(Elev), data=scores))
summary(nmds3.mod <- lm(MDS3 ~ Fish + scale(Elev), data=scores))

tab_model(nmds1.mod, nmds2.mod, nmds3.mod, 
          show.ci = F, show.se = T) # ,
        #  file = "./tables/nmds_axis_fit_mods_4Jun2025.doc")

nmds1_elev <- ggplot(scores) +
  geom_point(aes(x=Elev, y=MDS1)) +
  geom_smooth(aes(x=Elev, y=MDS1), color="black",method = "lm", se = FALSE) +
#  scale_color_manual(values=c("black") ) +
  labs(x="Elevation", y="NMDS1 Score", color=NULL) +
  geom_label(aes(x=3500, y=0.9), label = "\u03B2[Elev] = -0.49, \n R-sq = 0.66, p<0.0001", size = 4, family="Arial", label.size = 0) +
 # scale_color_manual(values = fishcolors) +
  theme_pubr() +
  theme(text = element_text(size=12, family = "Arial"), 
        axis.text = element_text(size=10), legend.position = "none")
#relevel(scores$Fish, "fish")

nmds2_fish <- ggplot(scores) + 
  geom_boxplot(aes(x=Fish, y=MDS2)) +
  ggpubr::geom_bracket( 
    tip.length = 0.02, # the downward "tips" of the bracket
    vjust = 0, # moves your text label (in this case, the p-value)
    xmin = 1, #starting point for the bracket
    xmax = 2, # ending point for the bracket
    y.position = 0.9, # vertical location of the bracket
    label.size = 4,# size of your bracket text
    label = paste0("\u03B2[Fishless] = -0.40, R-sq = 0.31 \n p<0.001"), family = "Arial") +
  labs(x="Lake Type", y="NMDS2 Score") +
  scale_y_continuous(limits = c(0,1.1)) +
  theme_pubr() +
  theme(text = element_text(size=12, family = "Arial"), 
        axis.text = element_text(size=10))

nmds.score.plot <- ggarrange(nmds1_elev, nmds2_fish, labels = c("(b)", "(c)"), font.label = list(size=12), vjust = 0)

all.NMDS <- ggarrange(NMDSplot, nmds.score.plot, ncol=1, heights = c(3,2), labels = c("(a)", ""), font.label = list(size=12))

 #ggsave(all.NMDS, filename = paste0("images/", "Figure3_", Sys.Date(), ".tiff"), 
  #      width = 18, height = 20, units = "cm", dpi = 600)


# Species ZINB ------------------------------------------------------------

allspp <- d3 %>% filter(location !="CENTER2R") %>% 
  group_by(year, basin, location, fish, date, jday, visit, loc.pt, Time, wind, species) %>%
  summarise(count = n()) %>%
  mutate(species = as.factor(species)) %>%
  left_join(lakes, by=c("location"="LAKENAME"))

allspp$fish <- factor(allspp$fish, levels=c("fishless", "fish"))
allspp$year <- factor(allspp$year, levels=c("2020", "2015", "2014"))

sppdets <- d3 %>% filter(location !="CENTER2R") %>% 
  group_by(year, basin, location, fish, date, jday, visit, loc.pt, time, wind, species) %>%
  summarise(count = n()) %>%
  group_by(species) %>% summarise(ndet=sum(count), nsite = n_distinct(location)/39) %>% arrange(desc(nsite)) %>%
  left_join(aoufull, by=c("species"="SPEC"))

keep <- sppdets %>% filter(ndet > 30) %>% arrange(species)

top9 <- allspp[allspp$species %in% keep$species,]

top9 %>% group_by(location, year) %>% summarise(nsurv = n_distinct(date)) #%>% View()

zerofill <- allspp %>% filter(species %in% keep$species) %>% 
  pivot_wider(names_from = species, values_from = count, values_fill = 0) %>%
  pivot_longer(GCRF:YRWA, names_to = "species", values_to = "count") %>% as.data.frame() %>%
  dplyr::select(year, basin, location, loc.pt, jday, Time, fish, ELEVATION, species, count)

zerofill$species <- factor(zerofill$species)

summary(allspp_pois <- glmer(count ~ fish*species + scale(ELEVATION) + (1|year) + (1|basin/location/loc.pt), data = zerofill, family = poisson))

overdisp_fun(allspp_pois) # indicates overdispersion 

# Compare model structures for zero-inflated component

summary(allspp_zip <- glmmTMB(count ~ fish*species + (1|basin/location/loc.pt),
                               zi = ~ 1, data = zerofill, family = poisson))

summary(allspp_zip2 <- glmmTMB(count ~ fish*species + scale(ELEVATION) + (1|basin/location/loc.pt),
                               zi = ~ scale(ELEVATION), data = zerofill, family = poisson))

summary(allspp_zip_se <- glmmTMB(count ~ fish*species + (1|basin/location/loc.pt),
                               zi = ~ species*scale(ELEVATION), data = zerofill, family = poisson)) # the best structure for the data, but not enough of it (singular fit)

summary(allspp_zip3 <- glmmTMB(count ~ fish*species + (1|basin/location/loc.pt),
                               zi = ~ scale(jday), data = zerofill, family = poisson) )

summary(allspp_zip4 <- glmmTMB(count ~ fish*species + (1|basin/location/loc.pt),
                               zi = ~ scale(Time), data = zerofill, family = poisson) )

anova(allspp_zip, allspp_zip2, allspp_zip3, allspp_zip4, allspp_zip_se) 

zerofill$species <- factor(zerofill$species)

plot(ggpredict(allspp_zip, terms=c("species", "fish"), type = "fe.zi")) # glance at marginal means preds

tab_model(allspp_zip, 
          transform = NULL,
          show.se = TRUE,
          show.df = T,
          show.fstat = T,
          dv.labels = c("model of count"),
    #      pred.labels = c("Intercept", "Fish"),
          string.se = "SE",
          string.ci = "95% CI",
          string.est = "Estimate",
          p.style = "numeric_stars",
          show.icc = F) #(),
      #    file = "./tables/ZIP_9spp_2025Jun4.doc") 

# retrieve marginal means +/- 95% CI for fish effect
zip.fx <- ggpredict(allspp_zip,terms = c("species", "fish"))
plot(predict_response(allspp_zip, terms = c("species", "fish"), type = "fe.zi", ci_level = 0.95)) +
 scale_color_manual(values=c(cbPalette[6], cbPalette[2])) +
  theme_pubr() 

modsumm <- summary(allspp_zip)

contrasts <- modsumm$coefficients$cond[11:18,]
contrasts <- as.data.frame(rbind(contrasts, modsumm$coefficients$cond[1,]) )# add AMPI
rownames(contrasts) <- c("CLNU", "DEJU", "DUFL", "GCRF", "HETH", "ROWR", "WCSP", "YRWA", "AMPI")
contrasts$species <- rownames(contrasts)

fx_contrasts <- left_join(as.data.frame(zip.fx), contrasts, by = c("x"="species"))

zip.plot <- ggplot(fx_contrasts) +
  geom_point(position=position_dodge(width=0.8), aes(x=reorder(x, `Pr(>|z|)`), y=predicted, color=group)) +
  geom_linerange(position=position_dodge(width=0.8),aes(x=reorder(x, `Pr(>|z|)`), ymin=conf.low, ymax=conf.high, color=group)) +
  ggpubr::geom_bracket(
    tip.length = 0.02, # the downard "tips" of the bracket
    vjust = 0, # moves your text label (in this case, the p-value)
    xmin = c(0.5, 1.5, 2.6), #starting point for the bracket
    xmax = c(1.4, 2.4, 3.4),# ending point for the bracket
    y.position = c(1.6, 3, 2.4), # vertical location of the bracket
    label.size = 3,# size of your bracket text
    label = c(paste0("p = 0.019"), paste0("p = 0.041"),paste0("p = 0.072")) ) +
  scale_color_manual(values=fishcolors) +
  #  geom_text(aes(x=1.5, y=10, label = paste("p = ", format(p_value, digits=2)) )) +
  theme_pubr() +
  labs(x=NULL, y="predicted per-point count", colour = "Lake Type") +
  theme(text = element_text(size=12, family = "Arial"), 
        axis.text = element_text(size=10))
zip.plot
#ggsave(zip.plot, filename = paste0("images/", "Figure4_", Sys.Date(), ".tiff"), 
#       width = 18, height = 12, units = "cm", dpi = 600)

# BACI --------------------------------------------------------------------

multraw <- d3 %>% filter(basin=="Center" | basin=="Amphitheater", location !="AMPHIT1B") %>% 
  group_by(year, basin, location, fish, date, jday, visit, point) %>% summarise(totabun = n(), totrich = n_distinct(species)) 

multraw %>% group_by(location, year) %>% summarise(min=min(date), max=max(date))

multraw$year <- factor(multraw$year)
multraw$location[multraw$location=="CENTER2R"] <- "CENTER2"
multraw$location <- factor(multraw$location)
multraw$basin <- factor(multraw$basin)
multraw$impact <- factor(if_else(multraw$location=="CENTER2", "removal", "control"))

multraw$impact2 <- factor(if_else(multraw$location=="CENTER2", "removal", if_else(multraw$location=="CENTER1", "inbasin_control", "outbasin_control")), levels=c("outbasin_control", "inbasin_control", "removal")) 
multraw$removal <- factor(ifelse(multraw$year != 2020, "pre", "post"), levels = c("pre", "post"), labels = c("Before ('14-'15)", "After ('20)"))
multraw$ba <- factor(if_else(multraw$removal=="Before ('14-'15)", "before", "after"), levels=c("before", "after"))
multraw$loc.pt <- paste(multraw$location, multraw$point, sep=".")
multraw$facet_labels <- ifelse(multraw$basin=="Amphitheater", "Control Basin", "Impact Basin")
multraw$fish[multraw$location=="CENTER2"] <- "fish removed"
multraw$fish <- factor(multraw$fish, levels=c("fishless", "fish", "fish removed"))

# Abundance
baci.pois.a <- glmer(totabun ~ ba*impact2 + (1|year) + (1|location/loc.pt),
                     family= poisson, data = multraw,
                     glmerControl(calc.derivs = F, optCtrl=list(maxfun=10^4)))
summary(baci.pois.a)
plot(ggpredict(baci.pois.a, terms=c("ba", "impact2")) )

# Richness
baci.pois.r <- glmer(totrich ~ ba*impact2 + (1|year) + (1|location/loc.pt),
                     family= poisson, data = multraw, 
                     glmerControl(calc.derivs = F, optCtrl=list(maxfun=10^4)))
summary(baci.pois.r)
plot(ggpredict(baci.pois.r, terms=c("ba", "impact2")) )

overdisp_fun(baci.pois.r) # good
overdisp_fun(baci.pois.a) # good

baci.model.tab <-
  tab_model(baci.pois.a, baci.pois.r,
            transform = NULL,
            show.se = FALSE,
            dv.labels = c("Abundance", "Richness"),
        #    pred.labels = c("Intercept", "Before-After [After]", "Control-Impact [Control Basin]", "Control-Impact [Impact]", "BA*CI [Control Basin]", "BA*CI [Impact]"),
         #   string.se = "SE",
            string.ci = "95% CI",
            string.est = "Estimate",
            p.style = "numeric",
            show.r2 = T,
            show.icc = F,
            show.aic = F,
            show.loglik = T,
            show.dev = T )#,
          #  file = "./tables/BACI_abun_richMods_2025Jun7.doc") 

# PLOTTING BACI 

mult <- d3 %>% filter(basin=="Center" | basin=="Amphitheater", location !="AMPHIT1B") %>% group_by(year, basin, location, fish, point, date, jday, visit) %>% summarise(totabun = n(), bypt_abun = n()/n_distinct(point), rich = n_distinct(species), bypt_rich = n_distinct(species)/n_distinct(point))
mult$removal <- factor(ifelse(mult$year == 2020, "post", "pre"), levels = c("pre", "post"), labels = c("Before \n('14-'15)", "After \n('20)"))
mult$location[mult$location=="CENTER2R"] <- "CENTER2"
mult$year <- factor(mult$year, levels=c(2020, 2015, 2014))
mult$loc.pt <- paste(mult$location, mult$point, sep=".")
mult_summary <- mult %>% 
  group_by(basin, location, fish, removal) %>% 
  summarise(n_pt = n_distinct(point),
            avg_abun=mean(totabun), 
            sd_abun=sd(totabun), 
            se_abun=sd_abun/sqrt(n_pt), 
            avg_rich=mean(rich), 
            sd_rich=sd(rich), 
            se_rich = sd_rich/sqrt(n_pt))

mult_summary$impact <- factor(ifelse(mult_summary$location=="CENTER2", "fish removed", "control lake"))
mult_summary$fish[mult_summary$location=="CENTER2"] <- "fish removed"
mult_summary$fish <- factor(mult_summary$fish, levels=c("fishless", "fish", "fish removed"))
mult_summary$facet_labels <- ifelse(mult_summary$basin=="Amphitheater", "Control Basin", "Impact Basin")


head(mult_summary)
mult$facet_labels <- ifelse(mult$basin=="Amphitheater", "Control Basin", "Impact Basin")
mult$fish[mult$location=="CENTER2"] <- "fish removed"
mult$fish <- factor(mult$fish, levels=c("fishless", "fish", "fish removed"))
mult$removal <- factor(mult$removal, levels=c("Before \n('14-'15)","After \n('20)"  ))
# ggplot(mult, aes(x=removal, y=totabun, color=fish)) +
#   geom_boxplot(position=position_dodge(width=1), aes(x=removal, y=totabun, color=fish)) +
#  # geom_line(position=position_dodge(width=1), aes(removal, totabun, group=location)) +
#   scale_color_manual(values = c(cbPalette[6], cbPalette[2], cbPalette[8]), name = "Lake Type") +
#   facet_grid(~facet_labels, scales="free") +
#   theme_pubr() +
#   labs(x=NULL, y="# of detected individuals", colour = "Lake Type") +
#   theme(text = element_text(size=12, family = "Times New Roman"),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# 
df <- mult %>% group_by(basin, location, fish, removal) %>%
  summarise(
  n_pt = n_distinct(point),
  y0 = min(bypt_rich),
  y25 = quantile(bypt_rich, 0.25),
  y50 = mean(bypt_rich),
  y75 = quantile(bypt_rich, 0.75),
  y100= max(bypt_rich)
) 

df$id <- paste(df$location, df$removal)
df$facet_labels <- ifelse(df$basin=="Amphitheater", "Control Basin", "Impact Basin")

ggplot(df) +
  geom_boxplot(aes(x = fish, ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100, color=removal),
    stat = "identity") +
  facet_wrap(~facet_labels, scales="free_x") +
  scale_color_manual(values = c(cbPalette[6], cbPalette[2], cbPalette[8]), name = "Lake Type") +
  theme(text = element_text(size=12, family = "Arial"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

baci.a.plot <- ggplot(mult_summary, aes(removal, avg_abun, color=fish)) +
  geom_crossbar(position=position_dodge(width=1), aes(ymin=avg_abun-se_abun, ymax=avg_abun+se_abun)) +
  geom_line(position=position_dodge(width=1), aes(removal, avg_abun, group=location)) +
  scale_color_manual(values = c(cbPalette[6], cbPalette[2], cbPalette[8]), name = "Lake Type") +
  facet_grid(~facet_labels, scales="free") +
  theme_pubr() +
  labs(x=NULL, y="no. of detected individuals", colour = "Lake Type") +
  theme(text = element_text(size=12, family = "Arial"), axis.text = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.background =element_rect(fill="white"))

#ggsave(plot=baci.a.plot, filename = "images/baci_amphit_center_abundance_2025Jun4.png", 
#       width = 6.5, height = 3, units = "in", dpi = 600)

ggpredict(baci.pois.r, terms=c("ba", "impact2"))

baci.r.plot <- ggplot(mult_summary, aes(removal, avg_rich, color=fish)) +
  geom_crossbar(position=position_dodge(width=1), aes(ymin=avg_rich-se_rich, ymax=avg_rich+se_rich)) +
  geom_line(position=position_dodge(width=1), aes(removal, avg_rich, group=location)) +
  theme_pubr() +
  scale_color_manual(values = c(cbPalette[6], cbPalette[2], cbPalette[8]), name = "Lake Type") +
  labs( x = NULL, y = "no. of detected species") +
  facet_grid(~facet_labels, scales="free") +
  theme(text = element_text(size=12, family = "Arial"), axis.text = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.background =element_rect(fill="white"))

#ggsave(plot=baci.r.plot, filename = "images/baci_amphit_center_richness_2025Jun4.png", 
#       width = 6.5, height = 3, units = "in", dpi = 600)

both.baci <- ggarrange(baci.a.plot, baci.r.plot, 
                       labels = c("(a)", "(b)"), 
                       ncol = 2, nrow = 1, font.label = list(size=12),
                       common.legend = T) 



# Shifts in species after fish removal ------------------------------------

center.raw <- d3 %>% filter(basin=="Center", jday > 158) %>% group_by(year, location, date, jday, point, visit, species) %>% summarise(count=n())
center.raw$location[center.raw$location=="CENTER2R"] <- "CENTER2"
center.raw %>% filter(location=="CENTER2", year != 2020) %>% group_by(year) %>% summarise(surveys = n_distinct(date)) # 5 unique survey dates

# LAKE LEVEL COUNTS FOR GOLDENBEAR
gb_pre <- center.raw %>% filter(location=="CENTER2", year != 2020) %>%
  group_by(species,date) %>% summarise(lakeabun=sum(count)) %>%
  group_by(species) %>% summarise(n_surv=n_distinct(date), meanabun=mean(lakeabun), med = median(lakeabun), sd_abun=sd(lakeabun), se_abun=sd_abun/sqrt(n_surv))

# with zeros
gb_pre <- center.raw %>% filter(location=="CENTER2", year != 2020) %>%
  group_by(species,date) %>% summarise(lakeabun=sum(count)) %>%
  spread(species, lakeabun, fill=0) %>%
  pivot_longer(AMPI:YRWA, values_to="count", names_to="species") %>%
  group_by(species) %>% summarise(nsurv = n_distinct(date), meanabun=mean(count), med = median(count), sd_abun=sd(count), se_abun=sd_abun/sqrt(nsurv))

gb_post <- center.raw %>% filter(location=="CENTER2", year == 2020) %>% group_by(species) %>% summarise(abun=sum(count))

gb_prepost <- full_join(gb_pre, gb_post, by="species") %>% mutate(postabun = abun)
gb_prepost[is.na(gb_prepost)] <- 0
gb_prepost$abun <- NULL
#write_csv(gb_prepost, "data/gb_prepost_byspp.csv")

gb_prepost$isup <- ifelse(gb_prepost$postabun > (gb_prepost$meanabun+gb_prepost$se_abun), "up","no")
gb_prepost$isup <- ifelse(gb_prepost$postabun >= gb_prepost$meanabun-gb_prepost$se_abun & gb_prepost$postabun <= (gb_prepost$meanabun+gb_prepost$se_abun), "in range", gb_prepost$isup)
gb_prepost$isup <- ifelse(gb_prepost$postabun < (gb_prepost$meanabun-gb_prepost$se_abun), "down",gb_prepost$isup)
gb_prepost$isup <- ifelse(gb_prepost$se_abun==0, "nd", gb_prepost$isup)


gb_prepost$isup[gb_prepost$species=="AMDI" | gb_prepost$species=="MOBL"] <- "new species"

gb_prepost$isup <- factor(gb_prepost$isup, levels=c("in range", "up", "new species"))

prepostspp <- ggplot(gb_prepost) + 
  geom_crossbar(aes(x=species, y=meanabun, ymin=meanabun-se_abun, ymax=meanabun+se_abun)) +
  geom_point(aes(x=species, y=postabun, color=isup), size=3) + 
  geom_linerange(aes(x=species, ymin=meanabun, ymax=postabun, color=isup)) +
  theme_classic() + 
  scale_color_manual(values=c(cbPalette[1], cbPalette[3], cbPalette[4])) +
  labs(y="count of individuals",x="Species Code", color="change") +
  scale_y_continuous(breaks=seq(0, 10, 1)) +
  theme(legend.position = "top",
        text = element_text(size=12, family = "Arial"), axis.text = element_text(size=10),
        axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


baciplots <- ggarrange(both.baci, prepostspp,
                       labels = c("", "(c)"), font.label = list(size=12),
                       ncol = 1, nrow = 2, 
                       common.legend = F) 
# ggsave(plot = baciplots, filename = paste0("images/", "Figure6_", Sys.Date(), ".tiff"),
#       width = 18, height = 22, units = "cm", dpi = 600)


# Mayflies ----------------------------------------------------------------

mayfly$roundw <- ifelse(mayfly$round==1, "June", ifelse(mayfly$round==2, "July", "August"))
mayfly$roundw <- factor(mayfly$roundw, levels = c("June", "July", "August"))
mayfly$fish <- factor(ifelse(mayfly$fish=="fish-containing", "stocked", "fishless"))

summary(zinb.mayfly <- glmmTMB(mayfly ~ fish + roundw + (1|basin/lake),
                               zi = ~ fish, offset = daysin, data=mayfly, family = poisson))

mayfly %>% left_join(lakes, by=c("lake"="LAKENAME")) %>% group_by(basin, LAKE, fish, roundw) %>% summarise(n = sum(mayfly)) %>%
  pivot_wider(names_from = roundw, values_from = n) %>% write.csv("./tables/mayfly_count_table.csv")


mayflies_round <- ggplot(mayfly) + 
  geom_jitter(position=position_jitterdodge(jitter.width = 0.4), aes(x=roundw, y=mayfly, color=fish, fill=fish),alpha=0.6) +
  geom_boxplot(aes(x=roundw, y=mayfly, fill=fish, color=fish), alpha=0.6) +
  scale_fill_manual(values = c(cbPalette[6], cbPalette[2])) +
  scale_color_manual(values = c(cbPalette[6], cbPalette[2])) +
  labs(title = NULL, 
       x = "Sampling Round",
       y = "Mayfly Count",
       color ="Lake Type",
       fill="Lake Type") + 
  theme_bw() +
  theme(legend.position="top", 
        axis.text = element_text(size=10),
        text = element_text(color="black", size=12),
        strip.text = element_text(color = "black", size = 12),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Additional Sticky trap data
bugs <- read_csv("data/insects/Results.csv") %>% 
  separate(Label, c("lake", "dir", "datein", "extra"), sep = "_") %>% 
  mutate(dateout=ymd(substr(extra, 1,8)), datein=ymd(datein), daysin=dateout-datein, sample = `...1`) %>%
  dplyr::select(sample, lake, dir, datein, dateout, daysin, Area) 
bugs$lake <- substr(bugs$lake, 1,7)
bugs$fish <- as.factor(substr(bugs$lake, 7,7))
bugs$basin <- substr(bugs$lake, 1,6)

# add a column for round and merge it with dataset
# 'bugs' is a csv where each row represents an individual on the imageJ scan. 
roundinfo <- bugs %>% group_by(lake, basin, fish, datein) %>% 
  summarise(totbugs = n()) %>%
  arrange(datein, .by_group=TRUE) %>% 
  mutate(round=as.factor(seq_along(datein))) %>% 
  dplyr::select(lake, datein, round)
bugs2 <- left_join(bugs, roundinfo) 
bugs2$datein <- ymd(bugs2$datein)
bugs2$daysin <- as.numeric(bugs2$daysin)
bugs2$roundw <- as.factor(ifelse(bugs2$round==1, "June", "July"), levels=c("June", "July") )
bugs2$Fish <- ifelse(bugs2$fish==1, "fishless", "stocked")
# summarise by individual plate
bugsplate <- bugs2 %>% 
  group_by(lake, basin, fish, round, dir, daysin) %>% 
  summarise(totbugs = n(), totarea = sum(Area), mean_bugsize = mean(Area)) %>%
  mutate(roundw = factor(ifelse(round==1,"June", "July"), levels=c("June", "July")),
         Fish = ifelse(fish==1, "fishless", "stocked"))

bugsizeplot <- ggplot(bugsplate) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3),aes(x=roundw, y=mean_bugsize, color=Fish, fill=Fish), alpha=0.6) +
  geom_boxplot(aes(x=roundw, y=mean_bugsize, fill=Fish, color=Fish), alpha=0.6) +
  scale_fill_manual(values=c(cbPalette[6], cbPalette[2]))  +
  scale_color_manual(values=c(cbPalette[6], cbPalette[2]))  +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y=expression(paste("mean insect size (mm", ""^2, ")")), x=NULL) +
  theme(legend.position="none", 
        text = element_text(color="black", size=12, family="Arial"),
        axis.text = element_text(size=10),
        strip.text = element_text(color = "black", size = 10),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

areaplot <- ggplot(bugsplate) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3), aes(x=roundw, y=totarea, color=Fish, fill=Fish), alpha=0.6) +
  geom_boxplot(aes(x=roundw, y=totarea, color=Fish, fill=Fish), alpha=0.6) +
  scale_fill_manual(values = c(cbPalette[6], cbPalette[2])) +
  scale_color_manual(values = c(cbPalette[6], cbPalette[2])) +
  theme_bw() +
  labs(y = expression(paste("trap area covered by insects (mm", ""^2, ")")), x=NULL) +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(color="black", size=12, family="Arial"),
        axis.text = element_text(size=10),
        strip.text = element_text(color = "black", size = 10),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

bugcountplot <- ggplot(bugsplate) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3), aes(x=roundw, y=totbugs, color=Fish, fill=Fish), alpha=0.6) +
  geom_boxplot(aes(x=roundw, y=totbugs, fill=Fish, color=Fish), alpha=0.6) +
  scale_fill_manual(values=c(cbPalette[6], cbPalette[2]))  +
  scale_color_manual(values=c(cbPalette[6], cbPalette[2]))  +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y="no. of insects per trap", x=NULL) +
  theme_bw() +
  theme(legend.position="bottom", 
        text = element_text(color="black", size=12, family="Arial"),
        axis.text = element_text(size=10),
        strip.text = element_text(color = "black", size = 10),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# # percent area covered by trap
# 
# plates <- read_csv("data/insects/Summary1.csv") %>% 
#   separate(Slice, c("Lake", "dir", "datein", "extra"), sep = "_") %>% 
#   mutate(dateout=ymd(substr(extra, 1,8)), datein=ymd(datein),daysin=dateout-datein) %>%
#   dplyr::rename(count=Count, areacov = `Total Area`, avgsize = `Average Size`, pctarea = `%Area`) %>%
#   dplyr::select(Lake, dir, datein, daysin, count, areacov, pctarea, avgsize)
# 
# plates$Lake <- substr(plates$Lake, 1,7)
# plates$fish <- as.factor(substr(plates$Lake, 7,7))
# 
# plates2 <- plates %>% left_join(roundinfo)
# 
# plates2$countperday <- plates2$count / as.numeric(plates2$daysin)
# 
# perc_area <- ggplot(plates2) +
#   geom_jitter(position=position_jitterdodge(jitter.width = 0.3), aes(x=round, y=pctarea, color=fish, fill=fish), alpha=0.6) +
#   geom_boxplot(aes(x=round, y=pctarea, color=fish, fill=fish), alpha=0.6) +
#   scale_fill_manual(values = c(cbPalette[6], cbPalette[2])) +
#   scale_color_manual(values = c(cbPalette[6], cbPalette[2])) +
#   theme_bw() +
#   labs(y = "% trap area covered", x=NULL) +
#   theme_bw() +
#   theme(legend.position="none", 
#         text = element_text(color="black", size=12, family="Arial"),
#         axis.text = element_text(size=10),
#         strip.text = element_text(color = "black", size = 10),
#         strip.background = element_rect(color="black", fill="white"),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())

bugs2$sampleID <- paste(bugs2$lake, bugs2$dir, bugs2$datein, sep="_")

sizedeets <- bugs2 %>% group_by(Fish, round) %>% summarise(n_plates = n_distinct(sampleID), mean_size = mean(Area), se_size = sd(Area)/sqrt(n_plates))


imagej_plots <- ggarrange(bugcountplot, areaplot, bugsizeplot, common.legend = T,
                          ncol=3, legend = "none")

allbugs <- ggarrange(mayflies_round, imagej_plots, nrow=2, heights = c(1.5,1))

# ggsave(plot = allbugs, filename = "images/all_insect_data_plots.png",
#        width = 6.5, height = 6, units = "in")

# Tables (Lakes, Species) -------------------------------------------------

head(lakes)
knappdb <- read_csv("data/knappdb_survey.csv")
head(knappdb)
laketab <- dat %>% group_by(year, basin, location, fish, date, visit) %>% summarise(npt= n_distinct(point))%>% spread(year, date) %>% left_join(lakes, by=c("location"="LAKENAME")) %>%
  dplyr::select(basin, LAKE, fish, NPOINTS, visit, `2014`, `2015`, `2020`)
write_csv(laketab, "./tables/laketables_draft.csv")

# species table
sppdets <- d3 %>% 
  mutate(locpoint = paste(location, point, sep = ".")) %>%
  group_by(year, basin, location, fish, date, jday, visit, locpoint, time, wind, species) %>%
  summarise(count = n()) %>%
  group_by(species) %>% summarise(ndet=sum(count), nsite = n_distinct(location)/39, fishcount = n_distinct(location[fish=="fish"]), fishlesscount = n_distinct(location[fish=="fishless"])) %>% arrange(desc(nsite)) %>%
  left_join(aoufull, by=c("species"="SPEC"))

sppdets$in_model <- ifelse(sppdets$species %in% unique(top9$species), "Y", "N")
arrange(sppdets, desc(nsite))
write_csv(sppdets, "./tables/specieslist_detailed.csv")

# lake figure
forplot <- lakes %>% dplyr::select(LAKE, FISH, ELEVATION, DEPTH, PERIMETER, NPOINTS) %>% 
  #mutate(`AREA km^2` = AREA/1e6) %>% 
  pivot_longer(cols = ELEVATION:NPOINTS)
forplot$FISH <- as.factor(ifelse(forplot$FISH=="1", "fishless (n=22)", "stocked (n=17)"))

forplot$facet_label[forplot$name=="ELEVATION"] <- "ELEVATION (m)"
forplot$facet_label[forplot$name=="DEPTH"] <- "MAXIMUM DEPTH (m)"
forplot$facet_label[forplot$name=="PERIMETER"] <- "PERIMETER (m)"
forplot$facet_label[forplot$name=="NPOINTS"] <- "NUMBER of BIRD SURVEY POINTS"


env.covar <- ggplot(forplot, aes(x=value, fill=FISH)) +
  geom_density(color="black", alpha=0.6, position = "identity") +
  scale_fill_manual(values=c(cbPalette[6], cbPalette[2])) +
  theme_pubclean() +
  facet_wrap(~facet_label, scales = "free") +
  labs(fill = "presence/absence of introduced fish") +
  theme(text = element_text(size=10, family = "Helvetica"), axis.title = element_blank())


#ggsave(plot=env.covar, filename = "images/environ_covars_by_lake_20250602.png", width = 6.5, height = 5, units = "in", dpi=600)

lakes %>% group_by(FISH) %>% summarise(tot_perim = sum(PERIMETER), tot_survs = sum(NPOINTS), scaled = tot_perim/tot_survs)
