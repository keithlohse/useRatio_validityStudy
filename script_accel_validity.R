library(tidyverse); library(ggpubr); library(ggpmisc); library(lme4);
library(car); library(lmerTest); library(GGally); library(viridis); library(patchwork)


# Set Colorblind friendly pallette 
cbPalette <- c("#999999", "#56B4E9","#E69F00", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#999933", "#882255", "#661100", "#6699CC")

# Set Colorblind friendly pallette 
cbGradient <- c("#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda","#8c96b9",
                "#8c96c6", "#8c6bb1", "#88419d", "#810f7c", "#4d004b")

setwd("C:/Users/kelop/Box/LangLab/DataHarmonization/useRatio_validityStudy")
list.files()

UE_DATA <- read.csv("./UpperLimbAccelerometry_restricted.csv", header=TRUE, stringsAsFactors = TRUE)
LE_DATA <- read.csv("./LowerLimbAccelerometry_restricted.csv", header=TRUE, stringsAsFactors = TRUE)
CLINICAL <- read.csv("./DemographicClinicalData_restricted.csv", header=TRUE, stringsAsFactors = TRUE)

colnames(CLINICAL)
colnames(UE_DATA)
colnames(LE_DATA)

STROKE <- CLINICAL %>%
  filter(Condition==1) %>%
  mutate(SubIDName=factor(SubIDName))

LE_DATA <- LE_DATA %>% 
  pivot_wider(values_from = AverageStepsPerDay, 
              names_from = SensorLocation, names_prefix = "Steps_")

STROKE <- merge(x=STROKE, y=UE_DATA, by=c("SubIDName", "TimePoint"))
STROKE <- merge(x=STROKE, y=LE_DATA, by=c("SubIDName", "TimePoint"),
                all = TRUE)


colnames(STROKE)
write.csv(data.frame(xtabs(~SubIDName+TimePoint, data=STROKE)) %>%
  pivot_wider(names_from = TimePoint, values_from = Freq), 
  "./outputs/data_at_each_time.csv")

summary(as.factor(STROKE$AffectedSide))

# Note that participants from "PMC8442937" all had first strokes, however,
# because this was an inclusion criterion, it was not captured in a REDCap field
# below, we change all of the NAs for these participants to 1's indicating their first stroke
str_sub(STROKE$SubIDName, 1, 10)

# Note that we filter out anyone who potentially had bilateral deficits below, but 
# no one actually does! Everyone in the stroke cohort either had a confirmed 
# dominant or nondominant side deficit.
STROKE %>% filter(AffectedSide==3) %>% group_by(SubIDName) %>% slice(1)
STROKE %>% filter(AffectedSide==4) %>% group_by(SubIDName) %>% slice(1)

# recode arms as preferred and non-preferred
STROKE <- STROKE %>%
  filter(AffectedSide!=3) %>% # exclude both sides affected
  filter(AffectedSide!=4) %>% # exclude neither side affected
  mutate(
    NumberofStrokes = ifelse(str_sub(SubIDName, 1, 10)=="PMC8442937", 
                             1, NumberofStrokes), # set all participants from PMC 8442937 to 1
    non_time = ifelse(AffectedSide == "2", l_time, r_time), #2 = right side affected
         par_time = ifelse(AffectedSide == "2", r_time,l_time),
         non_only_time = ifelse(AffectedSide == "2", l_only_time, r_only_time),
         par_only_time = ifelse(AffectedSide == "2", r_only_time, l_only_time),
         non_magnitude = ifelse(AffectedSide == "2", l_magnitude, r_magnitude),
         par_magnitude = ifelse(AffectedSide == "2", r_magnitude, l_magnitude),
         non_magnitude_sd = ifelse(AffectedSide == "2", l_magnitude_sd,r_magnitude_sd),
         par_magnitude_sd = ifelse(AffectedSide == "2", r_magnitude_sd,l_magnitude_sd),
         non_peak_magnitude = ifelse(AffectedSide == "2", l_peak_magnitude,r_peak_magnitude),
         par_peak_magnitude = ifelse(AffectedSide == "2", r_peak_magnitude,l_peak_magnitude),
         non_entropy = ifelse(AffectedSide == "2", l_entropy,r_entropy),
         par_entropy = ifelse(AffectedSide == "2", r_entropy,l_entropy),
         non_jerk_ave = ifelse(AffectedSide == "2", l_jerk_avg,r_jerk_avg),
         par_jerk_ave = ifelse(AffectedSide == "2", r_jerk_avg,l_jerk_avg),
         non_mean_freq = ifelse(AffectedSide == "2", l_mean_freq,r_mean_freq),
         par_mean_freq = ifelse(AffectedSide == "2", r_mean_freq,l_mean_freq),
         non_sd_freq = ifelse(AffectedSide == "2", l_sd_freq,r_sd_freq),
         par_sd_freq = ifelse(AffectedSide == "2", r_sd_freq,l_sd_freq),
  )

colnames(STROKE)
head(STROKE)


# Create Chronicity Variables -------------------------------------------------------
# correlations between variables over time
colnames(STROKE)
summary(STROKE$Chronicity)
# chronicity is in years, multiply by 52.1 to get approximate weeks
summary(STROKE$TimePoint)

LONG <- STROKE %>% select(SubIDName, Chronicity, TimePoint, StrokeType, 
                          AffARATTotal, UEFuglMeyer, 
                        use_ratio, par_time, non_time) %>%
  group_by(SubIDName, TimePoint) %>% # regroup by subject and then sort by ascending times
  slice(1) %>% # remove repeated rows from accel data to take first row at each time for each person
  group_by(SubIDName) %>% 
  mutate(enrollChronicity = Chronicity[1]*52.1) %>% # convert enrollment chronicity to weeks
  mutate(# convert years to weeks then add together:
    SubIDName = factor(SubIDName),
    WeeksPostStroke = enrollChronicity+TimePoint,
    enrollCat = factor(cut(enrollChronicity, breaks=c(0, 6, 12, 24, 36, 52, 9000),
                       labels = c("Week0", "Week6", "Week12", "Week24", "Week36", "MoreThan52"))),
    enrollCat = fct_relevel(enrollCat, "Week0", "Week6", "Week12", "Week24", "Week36", "MoreThan52"),
    weeksCat = factor(cut(WeeksPostStroke, breaks=c(0, 6, 12, 24, 36, 52, 9000),
                      labels = c("Week0", "Week6", "Week12", "Week24", "Week36", "MoreThan52"))),
    weeksCat = fct_relevel(weeksCat, "Week0", "Week6", "Week12", "Week24", "Week36", "MoreThan52")
    ) %>%
  group_by(SubIDName, weeksCat) %>% arrange(TimePoint) %>% 
  relocate(SubIDName, enrollChronicity, TimePoint, WeeksPostStroke, enrollCat, weeksCat) %>%
  ungroup()

# Availability of Data at Different Time Points ---------------------------------
library(ggridges)
colnames(LONG)
LONG %>% select(SubIDName, weeksCat, WeeksPostStroke, AffARATTotal,
                      UEFuglMeyer, use_ratio, par_time, 
                      non_time) %>%
  group_by(weeksCat) %>%
  pivot_longer(cols=AffARATTotal:non_time, 
               names_to = "variable",
               values_to = "value") %>%
  mutate(present = is.na(value)==FALSE,
         variable = factor(variable))


ggplot(LONG %>% select(SubIDName, weeksCat, WeeksPostStroke, AffARATTotal,
                             UEFuglMeyer, use_ratio, par_time, 
                             non_time) %>%
         group_by(weeksCat) %>%
         pivot_longer(cols=AffARATTotal:non_time, 
                      names_to = "variable",
                      values_to = "value") %>%
         mutate(present = is.na(value)==FALSE,
                variable = factor(variable)) %>%
         na.omit(), 
       aes(x = WeeksPostStroke, 
           y = variable)) +
  stat_density_ridges(aes(fill=variable), 
    alpha=0.5, jittered_points = TRUE,
    position = position_points_jitter(width = 0.00, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 0.5) +
  scale_x_log10(name="Weeks Post Stroke", breaks=c(1,6,12,24,36,52,100,1000))+
  scale_y_discrete(name = "Outcome")+
  scale_fill_manual(values=cbPalette)+
  theme_ridges() + 
  theme(axis.text=element_text(size=10, color="black"), 
        legend.text=element_text(size=10, color="black"),
        legend.title=element_text(size=10, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/data_overTime.jpeg",
  plot = last_plot(),
  width = 5,
  height = 3,
  units = "in",
  dpi = 300
)






# Part 1: Cross sectional/between subject correlations -------------------------
FIRST_DATA <- LONG %>% group_by(SubIDName, weeksCat) %>% arrange(TimePoint) %>% 
  slice(1) %>%
  filter(is.na(weeksCat)==FALSE) %>%
  ungroup()

xtabs(~SubIDName, data=FIRST_DATA) 
xtabs(~SubIDName+enrollCat, data=FIRST_DATA) # number of observations by when enrolled
xtabs(~SubIDName+weeksCat, data=FIRST_DATA) # there should only be one observations by weekCat

FIRST_DATA %>% group_by(weeksCat) %>%
  summarise(complete_UR = sum(is.na(use_ratio)==FALSE),
            complete_ARAT = sum(is.na(AffARATTotal)==FALSE),
            complete_FM = sum(is.na(UEFuglMeyer)==FALSE))
  


# ARAT -------------------------------------------------------------------------
# Spaghetti plot showing data for all participants ---
ggplot(data=LONG %>% filter(is.na(weeksCat)==FALSE), 
       aes(x=WeeksPostStroke, y=AffARATTotal))+
  geom_rect(xmin=log(1.4, 10), xmax=log(4, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(6, 10), xmax=log(7, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(12, 10), xmax=log(14, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(24, 10), xmax=log(27, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(36, 10), xmax=log(40, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(52, 10), xmax=log(58, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_line(aes(group=SubIDName), col="black", lty=1, alpha=0.3)+
  geom_point(aes(group=SubIDName, fill=weeksCat), col="black", shape=21, alpha=0.3)+
  # stat_smooth(aes(group=SubIDName, col=SubIDName), method = "lm", formula=y~x,
  #             alpha=0.5, se=FALSE, lwd=0.5)+
  scale_x_log10(name = "Weeks Post-Stroke (log[10] scale)") +
  scale_y_continuous(name = "ARAT (Affected Arm)") +
  scale_fill_manual(values=cbPalette)+
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+labs(fill="Time Window")+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=10, color="black"),
        legend.title=element_text(size=10, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/ARAT_overTime_allSubjects.jpeg",
  plot = last_plot(),
  width = 5,
  height = 3,
  units = "in",
  dpi = 300
)

# Spaghetti plot of Use Ratio for all subjects --------------------
ggplot(data=LONG %>% filter(is.na(weeksCat)==FALSE), 
       aes(x=WeeksPostStroke, y=use_ratio))+
  geom_rect(xmin=log(1.4, 10), xmax=log(4, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(6, 10), xmax=log(7, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(12, 10), xmax=log(14, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(24, 10), xmax=log(27, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(36, 10), xmax=log(40, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(52, 10), xmax=log(58, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_line(aes(group=SubIDName), col="black", lty=1, alpha=0.3)+
  geom_point(aes(group=SubIDName, fill=weeksCat), col="black", shape=21, alpha=0.3)+
  # stat_smooth(aes(group=SubIDName, col=SubIDName), method = "lm", formula=y~x,
  #             alpha=0.5, se=FALSE, lwd=0.5)+
  scale_x_log10(name = "Weeks Post-Stroke (log[10] scale)") +
  scale_y_continuous(name = "Use Ratio") +
  scale_fill_manual(values=cbPalette)+
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+labs(fill="Time Window")+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=10, color="black"),
        legend.title=element_text(size=10, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/useRatio_overTime_allSubjects.jpeg",
  plot = last_plot(),
  width = 5,
  height = 3,
  units = "in",
  dpi = 300
)




# Spaghetti Plot of FM UE for all subjects ---------------------
colnames(LONG)
ggplot(data=LONG %>% filter(is.na(weeksCat)==FALSE), 
       aes(x=WeeksPostStroke, y=UEFuglMeyer))+
  geom_rect(xmin=log(1.4, 10), xmax=log(4, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(6, 10), xmax=log(7, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(12, 10), xmax=log(14, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(24, 10), xmax=log(27, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(36, 10), xmax=log(40, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_rect(xmin=log(52, 10), xmax=log(58, 10), ymin=-Inf, ymax=Inf, 
            col="grey60", fill="grey90")+
  geom_line(aes(group=SubIDName), col="black", lty=1, alpha=0.3)+
  geom_point(aes(group=SubIDName, fill=weeksCat), col="black", shape=21, alpha=0.3)+
  # stat_smooth(aes(group=SubIDName, col=SubIDName), method = "lm", formula=y~x,
  #             alpha=0.5, se=FALSE, lwd=0.5)+
  scale_x_log10(name = "Weeks Post-Stroke (log[10] scale)") +
  scale_y_continuous(name = "Fugl-Meyer UE") +
  scale_fill_manual(values=cbPalette)+
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+labs(fill="Time Window")+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=10, color="black"),
        legend.title=element_text(size=10, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/fuglyMeyer_overTime_allSubjects.jpeg",
  plot = last_plot(),
  width = 5,
  height = 3,
  units = "in",
  dpi = 300
)


# Paretic to Non-Paretic Time ----
for (t in c("Week0", "Week6", "Week12", "Week24", "Week36", "MoreThan52")) {
  print(t)
  print(cor.test(y=FIRST_DATA[FIRST_DATA$weeksCat==eval(t),]$par_time, 
           x=FIRST_DATA[FIRST_DATA$weeksCat==eval(t),]$non_time,
           method="pearson",
           conf.level = 0.95,
           use="complete"))
}



ggplot(data=FIRST_DATA, aes(y=par_time, x=non_time))+
  geom_point(col="black", shape=21)+
  #stat_smooth(method="lm", se=TRUE)+
  stat_poly_line() +
  stat_poly_eq() +
  geom_abline(intercept = 0, slope=1, col="black", lwd=0.5)+
  scale_x_continuous(name="Non-Paretic Time (h)", breaks=c(seq(0,12,2)), limits=c(0,12))+
  scale_y_continuous(name="Paretic Time (h)", breaks=c(seq(0,12,2)), limits=c(0,12)) +
  facet_wrap(~weeksCat, ncol=1)+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")

ggsave(
  filename="./outputs/Paretic_to_nonParetic_overTime.jpeg",
  plot = last_plot(),
  width = 4,
  height = 18,
  units = "in",
  dpi = 300
)


# ARAT by Use Ratio ----
for (t in c("Week0", "Week6", "Week12", "Week24", "Week36", "MoreThan52")) {
  print(t)
  print(cor.test(y=FIRST_DATA[FIRST_DATA$weeksCat==eval(t),]$AffARATTotal, 
                 x=FIRST_DATA[FIRST_DATA$weeksCat==eval(t),]$use_ratio,
                 method="pearson",
                 conf.level = 0.95,
                 use="complete"))
}

ggplot(data=FIRST_DATA, aes(y=AffARATTotal, x=use_ratio))+
  geom_point(shape=21)+
  #stat_smooth(aes(col=as.factor(StrokeType)), method="lm", se=TRUE)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_color_manual(values=cbPalette)+
  scale_x_continuous(name="Use Ratio (paretic/non-paretic)")+
  scale_y_continuous(name = "Affected Side ARAT", breaks=c(seq(0,60,10)),
                     limits=c(0,60)) +
  facet_wrap(~weeksCat, ncol=1)+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")

ggsave(
  filename="./outputs/useRatio_ARAT_overTime.jpeg",
  plot = last_plot(),
  width = 4,
  height = 18,
  units = "in",
  dpi = 300
)



# ARAT by Paretic Time ----
ggplot(data=FIRST_DATA, aes(y=AffARATTotal, x=par_time))+
  geom_point(col="black", shape=21)+
  #stat_smooth(method="lm", se=TRUE)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name="Paretic Time (h)", breaks=c(seq(0,12,2)))+
  scale_y_continuous(name = "Affected Side ARAT", breaks=c(seq(0,60,10)),
                     limits=c(0,60)) +
  facet_wrap(~weeksCat, ncol=1)+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")

ggsave(
  filename="./outputs/paretic_ARAT_overTime.jpeg",
  plot = last_plot(),
  width = 4,
  height = 18,
  units = "in",
  dpi = 300
)



summary(FIRST_DATA$par_time)
summary(FIRST_DATA$non_time)

# ARAT by Non-Paretic Time ----
ggplot(data=FIRST_DATA, aes(y=AffARATTotal, x=non_time))+
  geom_point(col="black", shape=21)+
  #stat_smooth(method="lm", se=TRUE)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name="Non-Paretic Time (h)", breaks=c(seq(0,12,2)))+
  scale_y_continuous(name = "Affected Side ARAT", breaks=c(seq(0,60,10)),
                     limits=c(0,60)) +
  facet_wrap(~weeksCat, ncol=1)+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")


ggsave(
  filename="./outputs/nonParetic_ARAT_overTime_WIDE.jpeg",
  plot = last_plot(),
  width = 4,
  height = 18,
  units = "in",
  dpi = 300
)








# FM UE by Use Ratio ----
for (t in c("Week0", "Week6", "Week12", "Week24", "Week36", "MoreThan52")) {
  print(t)
  print(cor.test(y=FIRST_DATA[FIRST_DATA$weeksCat==eval(t),]$UEFuglMeyer, 
                 x=FIRST_DATA[FIRST_DATA$weeksCat==eval(t),]$use_ratio,
                 method="pearson",
                 conf.level = 0.95,
                 use="complete"))
}


summary(FIRST_DATA$use_ratio)
ggplot(data=FIRST_DATA, aes(y=UEFuglMeyer, x=use_ratio))+
  geom_point(shape=21)+
  #stat_smooth(aes(col=as.factor(StrokeType)), method="lm", se=TRUE)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_color_manual(values=cbPalette)+
  scale_x_continuous(name="Use Ratio (paretic/non-paretic)")+
  scale_y_continuous(name = "Fugl-Meyer UE", breaks=c(seq(0,70,10)),
                     limits=c(0,70)) +
  facet_wrap(~weeksCat, ncol=1)+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")

ggsave(
  filename="./outputs/useRatio_fuglMeyer_overTime.jpeg",
  plot = last_plot(),
  width = 4,
  height = 18,
  units = "in",
  dpi = 300
)



# FM UE by Paretic Time ----
ggplot(data=FIRST_DATA, aes(y=UEFuglMeyer, x=par_time))+
  geom_point(col="black", shape=21)+
  #stat_smooth(method="lm", se=TRUE)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name="Paretic Time (h)", breaks=c(seq(0,12,2)))+
  scale_y_continuous(name = "Fugl-Meyer UE", breaks=c(seq(0,70,10)),
                     limits=c(0,70)) +
  facet_wrap(~weeksCat, ncol=1)+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")

ggsave(
  filename="./outputs/paretic_fulgMeyer_overTime.jpeg",
  plot = last_plot(),
  width = 4,
  height = 18,
  units = "in",
  dpi = 300
)




# FM UE by Non-Paretic Time ----
ggplot(data=FIRST_DATA, aes(y=UEFuglMeyer, x=non_time))+
  geom_point(col="black", shape=21)+
  #stat_smooth(method="lm", se=TRUE)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name="Non-Paretic Time (h)", breaks=c(seq(0,12,2)))+
  scale_y_continuous(name = "Fugl-Meyer UE", breaks=c(seq(0,70,10)),
                     limits=c(0,70)) +
  facet_wrap(~weeksCat, ncol=1)+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")

ggsave(
  filename="./outputs/nonParetic_fuglMeyer_overTime.jpeg",
  plot = last_plot(),
  width = 4,
  height = 18,
  units = "in",
  dpi = 300
)






# Part 1: Demographic statistics for the overall cohort --------------------------------
sublist <- unique(FIRST_DATA$SubIDName) # Make sure we are only pulling subject IDS from the cross sectional cohor
colnames(STROKE)

UNIQUE<-STROKE %>% group_by(SubIDName) %>% arrange(TimePoint) %>%
  mutate(Observations=n()) %>% ungroup() %>% arrange(SubIDName, TimePoint) %>%
  group_by(SubIDName) %>% slice(1) %>%
  mutate(AffectedSide=factor(AffectedSide),
         TherpyRehabLocTyp = factor(TherpyRehabLocTyp),
         HandPrefType=factor(HandPrefType),
         BirthSexAssignTyp=factor(BirthSexAssignTyp),
         EthnUSACat=factor(EthnUSACat),
         RaceUSACat=factor(RaceUSACat),
         StrokeType=factor(StrokeType),  
         StrokeLocation=factor(StrokeLocation),
         EduYrCt = factor(EduYrCt),
         LivingStatus=factor(LivingStatus),
         Chronicity=Chronicity*52.1) %>%
  select(SubIDName, AgeVal, BirthSexAssignTyp, RaceUSACat, EthnUSACat, EduYrCt,       
         StrokeType, StrokeLocation, AffectedSide, TherpyRehabLocTyp,
         Chronicity, NumberofStrokes, LivingStatus, 
         # Enroll Variables
         CESDTotalScore, MesulamTotalErrors, MesulamDifference, MoCATotal,
         AffARATTotal, UEFuglMeyer, non_time, par_time, use_ratio, Observations)

colnames(UNIQUE)

UNIQUE %>% summary()

summary(as.factor(UNIQUE$NumberofStrokes))






# Part 2: Longitudinal/within subject correlations -----------------------------
ACUTE <- LONG %>% filter(enrollCat=="Week0" | enrollCat=="Week6") %>%
  mutate(SubIDName = factor(SubIDName))

# Spaghetti Plots showing change in variables over time ------------------------
head(LONG)

# Change in ARAT
ggplot(data=ACUTE, 
       aes(x=WeeksPostStroke, y=AffARATTotal))+
  geom_line(aes(group=SubIDName), col="black", lty=1, alpha=0.3,
            position = position_dodge(width=0.4))+
  geom_point(aes(group=SubIDName), col="black", shape=21, alpha=0.3,
             position = position_dodge(width=0.4))+
  #stat_smooth(method = "loess", col="blue", alpha=0.5, se=TRUE, lwd=1,
  #            xseq = seq(0,26, length=80))+
  scale_x_continuous(name = "Weeks Post-Stroke") +
  scale_y_continuous(name = "ARAT (Affected Arm)") +
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/acute_ARAT_detailed.jpeg",
  plot = last_plot(),
  width = 5,
  height = 3,
  units = "in",
  dpi = 300
)



# Change in Use Ratio ----------------------------------------------------------
ggplot(data=ACUTE, 
       aes(x=WeeksPostStroke, y=use_ratio))+
  geom_line(aes(group=SubIDName), col="black", lty=1, alpha=0.3,
            position = position_dodge(width=0.4))+
  geom_point(aes(group=SubIDName), col="black", shape=21, alpha=0.3,
             position = position_dodge(width=0.4))+
  # stat_smooth(aes(group=SubIDName, col=SubIDName), method = "lm", formula=y~x,
  #              alpha=0.5, se=FALSE, lwd=0.5)+
  scale_x_continuous(name = "Weeks Post-Stroke") +
  scale_y_continuous(name = "Use Ratio") +
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/acute_useRatio_detailed.jpeg",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = "in",
  dpi = 300
)






# Change in ARAT for single subject ----
unique(ACUTE$SubIDName)
# PMC8442937_049
# PMC8442937_017
colnames(ACUTE)

ggplot(data=ACUTE %>% 
         filter(SubIDName=="PMC8442937_049"), 
       aes(x=WeeksPostStroke, y=AffARATTotal))+
  geom_line(aes(group=SubIDName), col="black", lty=1, 
            position = position_dodge(width=0.4))+
  geom_point(aes(group=SubIDName), col="black", shape=21,
             position = position_dodge(width=0.4))+
  #stat_smooth(method = "lm", formula=y~x+I(x^2)+I(x^3),
  #            alpha=0.5, se=FALSE, lwd=1)+
  scale_x_continuous(name = "Weeks Post-Stroke", limits=c(0,26)) +
  scale_y_continuous(name = "ARAT (Affected Arm)", limits=c(0,60),
                     breaks = c(seq(0,60,10))) +
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/change_in_ARAT_singleSubject_B.jpeg",
  plot = last_plot(),
  width = 3,
  height = 2.5,
  units = "in",
  dpi = 300
)


# Change in FM UE for a single subject ----
ggplot(data=ACUTE %>% 
         filter(SubIDName=="PMC8442937_049"), 
       aes(x=WeeksPostStroke, y=UEFuglMeyer))+
  geom_line(aes(group=SubIDName), col="black", lty=1, 
            position = position_dodge(width=0.4))+
  geom_point(aes(group=SubIDName), col="black", shape=21,
             position = position_dodge(width=0.4))+
  #stat_smooth(method = "lm", formula=y~x+I(x^2)+I(x^3),
  #            alpha=0.5, se=FALSE, lwd=1)+
  scale_x_continuous(name = "Weeks Post-Stroke", limits=c(0,26)) +
  scale_y_continuous(name = "Fugl-Meyer UE", limits=c(0,70),
                     breaks = c(seq(0,70,10))) +
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/change_in_fuglMeyer_singleSubject_B.jpeg",
  plot = last_plot(),
  width = 3,
  height = 2.5,
  units = "in",
  dpi = 300
)





# Change in use ratio for single subject ----
unique(ACUTE$SubIDName)
ggplot(data=ACUTE %>% 
         filter(SubIDName=="PMC8442937_017"), 
       aes(x=WeeksPostStroke, y=use_ratio))+
  geom_line(aes(group=SubIDName), col="black", lty=1, 
            position = position_dodge(width=0.4))+
  geom_point(aes(group=SubIDName), col="black", shape=21,
             position = position_dodge(width=0.4))+
  #stat_smooth(method = "lm", formula=y~x+I(x^2)+I(x^3),
  #            alpha=0.5, se=FALSE, lwd=1)+
  scale_x_continuous(name = "Weeks Post Stroke", limits=c(0,26)) +
  scale_y_continuous(name = "Use Ratio (par/non-par)", limits=c(0,1)) +
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/change_in_useRatio_singleSubject_A.jpeg",
  plot = last_plot(),
  width = 3,
  height = 2.5,
  units = "in",
  dpi = 300
)



# Knots for splines-------------------------------------------------------------
# ARAT: Write a for-loop to select best fitting spline ----
ACUTE[ACUTE$TimePoint==max(ACUTE$TimePoint),]
summary(ACUTE$TimePoint)

ARAT_SPLINES <- data.frame()
# note that the for loop works with integers, so we will work with integer 
# weeks from enrollment for the placement of our knots

# We know we want to place two knots, one generally around weeks 1 to 8
# and the other sometime after week 8
knot1 <- c(seq(1,10,1))
knot2 <- c(seq(11,20,1))

counter=0
for(i in knot1){
  print(i)
  
  # Create the spline
  ACUTE$spline1 <- ifelse(ACUTE$WeeksPostStroke<i,
                          0, # If time < i, put 0
                          ACUTE$WeeksPostStroke-i) # If time is >/=i, put {time - i}
  
  for(j in knot2){
    print(j)
    
    counter=counter+1
    
    # Create the spline
    ACUTE$spline2 <- ifelse(ACUTE$WeeksPostStroke<j,
                            0, # If time < i, put 0
                            ACUTE$WeeksPostStroke-j) # If time is >/=i, put {time - i}
    
    # Define the model
    spline_mod <-lmer(AffARATTotal~
                        # Fixed-effects
                        1+WeeksPostStroke+spline1+spline2+
                        # Random-effects
                        (1+WeeksPostStroke+spline1+spline2|SubIDName), 
                      data=ACUTE, 
                      REML=FALSE,
                      control=lmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=5e5)))
    
    
    # Store the AIC, other parameters can be added
    ARAT_SPLINES[counter, 1] <- i
    ARAT_SPLINES[counter, 2] <- j
    ARAT_SPLINES[counter, 3] <- AIC(spline_mod)
    
  }
}

ARAT_SPLINES
ARAT_SPLINES[ARAT_SPLINES$V3==min(ARAT_SPLINES$V3),]
ARAT_SPLINES %>% arrange(V3) %>% slice(c(1:10))
# Weeks 7 and 15 are the absolute minimum for ARAT,
# but that is not appreciably different from 6 and 14 or 6 and 13.

ggplot(ARAT_SPLINES, aes(x=V1, y=V2, fill = V3))+
  geom_tile(col="black") +
  scale_fill_viridis(discrete=FALSE, trans = 'reverse')+
  #scale_fill_gradient2(low="white", mid="dodgerblue", high="firebrick",
  #                     midpoint=3390)+
  scale_x_continuous(name = "First Knot Location", breaks = c(seq(0,10,1))) +
  scale_y_continuous(name = "Second Knot Location", breaks = c(seq(10,20,1))) +
  labs(fill="AIC")+
  ggtitle("ARAT Splines")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "right")

ggsave(
  filename="./outputs/ARAT_spline_fits.jpeg",
  plot = last_plot(),
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)



# Use Ratio: Write a for-loop to select best fitting spline ----
ACUTE[ACUTE$TimePoint==max(ACUTE$TimePoint),]
summary(ACUTE$TimePoint)

UR_SPLINES <- data.frame()
# note that the for loop works with integers, so we will work with integer 
# weeks from enrollment for the placement of our knots

# We know we want to place two knots, one generally around weeks 1 to 8
# and the other sometime after week 8
knot1 <- c(seq(1,10,1))
knot2 <- c(seq(11,20,1))

counter=0
for(i in knot1){
  print(i)
  
  # Create the spline
  ACUTE$spline1 <- ifelse(ACUTE$WeeksPostStroke<i,
                          0, # If time < i, put 0
                          ACUTE$WeeksPostStroke-i) # If time is >/=i, put {time - i}
  
  for(j in knot2){
    print(j)
    
    counter=counter+1
    
    # Create the spline
    ACUTE$spline2 <- ifelse(ACUTE$WeeksPostStroke<j,
                            0, # If time < i, put 0
                            ACUTE$WeeksPostStroke-j) # If time is >/=i, put {time - i}
    
    # Define the model
    spline_mod <-lmer(use_ratio~
                        # Fixed-effects
                        1+WeeksPostStroke+spline1+spline2+
                        # Random-effects
                        (1+WeeksPostStroke+spline1+spline2|SubIDName), 
                      data=ACUTE, 
                      REML=FALSE,
                      control=lmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=5e5)))
    
    
    # Store the AIC, other parameters can be added
    UR_SPLINES[counter, 1] <- i
    UR_SPLINES[counter, 2] <- j
    UR_SPLINES[counter, 3] <- AIC(spline_mod)
    
  }
}

UR_SPLINES
UR_SPLINES[UR_SPLINES$V3==min(UR_SPLINES$V3),]
UR_SPLINES %>% arrange(V3) %>% slice(c(1:10))
# Weeks 5 and 12 give the absolute best times for UR.
# but that is not appreciably different from 5 and 11 or 6 and 11.

ggplot(UR_SPLINES, aes(x=V1, y=V2, fill = V3))+
  geom_tile(col="black") +
  scale_fill_viridis(discrete=FALSE, trans = 'reverse')+
  #scale_fill_gradient2(low="white", mid="dodgerblue", high="firebrick",
  #                     midpoint=-570)+
  scale_x_continuous(name = "First Knot Location", breaks = c(seq(0,10,1))) +
  scale_y_continuous(name = "Second Knot Location", breaks = c(seq(10,20,1))) +
  ggtitle("Use-Ratio Splines")+
  labs(fill="AIC")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "right")


ggsave(
  filename="./outputs/UR_spline_fits.jpeg",
  plot = last_plot(),
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)




# Fugl Meyer UE: Write a for-loop to select best fitting spline ----
ACUTE[ACUTE$TimePoint==max(ACUTE$TimePoint),]
summary(ACUTE$TimePoint)

FM_SPLINES <- data.frame()
# note that the for loop works with integers, so we will work with integer 
# weeks from enrollment for the placement of our knots

# We know we want to place two knots, one generally around weeks 1 to 8
# and the other sometime after week 8
knot1 <- c(seq(1,10,1))
knot2 <- c(seq(11,20,1))

counter=0
for(i in knot1){
  print(i)
  
  # Create the spline
  ACUTE$spline1 <- ifelse(ACUTE$WeeksPostStroke<i,
                          0, # If time < i, put 0
                          ACUTE$WeeksPostStroke-i) # If time is >/=i, put {time - i}
  
  for(j in knot2){
    print(j)
    
    counter=counter+1
    
    # Create the spline
    ACUTE$spline2 <- ifelse(ACUTE$WeeksPostStroke<j,
                            0, # If time < i, put 0
                            ACUTE$WeeksPostStroke-j) # If time is >/=i, put {time - i}
    
    # Define the model
    spline_mod <-lmer(UEFuglMeyer~
                        # Fixed-effects
                        1+WeeksPostStroke+spline1+spline2+
                        # Random-effects
                        (1+WeeksPostStroke+spline1+spline2|SubIDName), 
                      data=ACUTE, 
                      REML=FALSE,
                      control=lmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=5e5)))
    
    
    # Store the AIC, other parameters can be added
    FM_SPLINES[counter, 1] <- i
    FM_SPLINES[counter, 2] <- j
    FM_SPLINES[counter, 3] <- AIC(spline_mod)
    
  }
}

FM_SPLINES
FM_SPLINES[FM_SPLINES$V3==min(FM_SPLINES$V3),]
FM_SPLINES %>% arrange(V3) %>% slice(c(1:10))
# Weeks 5 and 12 give the absolute best times for UR.
# but that is not appreciably different from 5 and 11 or 6 and 11.

ggplot(FM_SPLINES, aes(x=V1, y=V2, fill = V3))+
  geom_tile(col="black") +
  scale_fill_viridis(discrete=FALSE, trans = 'reverse')+
  #scale_fill_gradient2(low="white", mid="dodgerblue", high="firebrick",
  #                     midpoint=-570)+
  scale_x_continuous(name = "First Knot Location", breaks = c(seq(0,10,1))) +
  scale_y_continuous(name = "Second Knot Location", breaks = c(seq(10,20,1))) +
  ggtitle("Fugl-Meyer Splines")+
  labs(fill="AIC")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "right")


ggsave(
  filename="./outputs/FM_spline_fits.jpeg",
  plot = last_plot(),
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)




# Codes for different time variables
colnames(ACUTE)
ACUTE$orthTime.1 <- poly(ACUTE$WeeksPostStroke, 3)[,1]
ACUTE$orthTime.2 <- poly(ACUTE$WeeksPostStroke, 3)[,2]
ACUTE$orthTime.3 <- poly(ACUTE$WeeksPostStroke, 3)[,3]
ACUTE$invTime.2 <- (ACUTE$WeeksPostStroke+1)^-2
ACUTE$invTime.3 <- (ACUTE$WeeksPostStroke+1)^-3

plot(x=ACUTE$WeeksPostStroke, y=ACUTE$orthTime.1)
summary(ACUTE$WeeksPostStroke)

ACUTE$knot05 <- ifelse(ACUTE$WeeksPostStroke<5,
                      0, # If time < k, put 0
                      ACUTE$WeeksPostStroke-5) # If time is >/=k, put time

ACUTE$knot11 <- ifelse(ACUTE$WeeksPostStroke<11,
                       0, # If time < k, put 0
                       ACUTE$WeeksPostStroke-11) # If time is >/=k, put time

plot(x=ACUTE$WeeksPostStroke, y=ACUTE$knot05)
plot(x=ACUTE$WeeksPostStroke, y=ACUTE$knot11)


# Fitting Slopes for the ARAT --------------------------------------------------
colnames(ACUTE)
# Intercept model
int_mod_ARAT<-lmer(AffARATTotal~
                     # Fixed-effects
                     1+
                     # Random-effects
                     (1|SubIDName), 
                   data=ACUTE, 
                   REML=FALSE,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=5e5)))
summary(int_mod_ARAT)


# linear model
lin_mod_ARAT<-lmer(AffARATTotal~
                     # Fixed-effects
                     1+WeeksPostStroke+
                     # Random-effects
                     (1+WeeksPostStroke|SubIDName), 
                   data=ACUTE, 
                   REML=FALSE,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=5e5)))
summary(lin_mod_ARAT)

# quadratic model
sq_mod_ARAT<-lmer(AffARATTotal~
                    # Fixed-effects
                    1+WeeksPostStroke+I(WeeksPostStroke^2)+
                    # Random-effects
                    (1+WeeksPostStroke+I(WeeksPostStroke^2)|SubIDName), data=ACUTE, REML=FALSE,
                  control=lmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=5e5)))

summary(sq_mod_ARAT)

# cubic model
cub_mod_ARAT<-lmer(AffARATTotal~
                     # Fixed-effects
                     1+WeeksPostStroke+I(WeeksPostStroke^2)+I(WeeksPostStroke^3)+
                     # Random-effects
                     (1+WeeksPostStroke+I(WeeksPostStroke^2)|SubIDName), data=ACUTE, REML=FALSE,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=5e5)))
summary(cub_mod_ARAT)
anova(int_mod_ARAT, lin_mod_ARAT, sq_mod_ARAT, cub_mod_ARAT)

# inverse square root
invSq_ARAT <-lmer(AffARATTotal~
                    # Fixed-effects
                    1+invTime.2+
                    # Random-effects
                    (1+invTime.2|SubIDName), data=ACUTE, REML=FALSE,
                  control=lmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=5e5)))
summary(invSq_ARAT)

# inverse cube root
invCube_ARAT <-lmer(AffARATTotal~
                      # Fixed-effects
                      1+invTime.3+
                      # Random-effects
                      (1+invTime.3|SubIDName), data=ACUTE, REML=FALSE,
                    control=lmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=5e5)))
summary(invCube_ARAT)
anova(int_mod_ARAT, lin_mod_ARAT, sq_mod_ARAT, 
      cub_mod_ARAT, invSq_ARAT, invCube_ARAT)



# spline model
spline_mod_ARAT <-lmer(AffARATTotal~
                         # Fixed-effects
                         1+WeeksPostStroke+knot05+knot11+
                         # Random-effects
                         (1+WeeksPostStroke+knot05+knot11|SubIDName), data=ACUTE, REML=FALSE,
                       control=lmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=5e5)))
summary(spline_mod_ARAT)

anova(int_mod_ARAT, lin_mod_ARAT, sq_mod_ARAT, 
      cub_mod_ARAT, invSq_ARAT, invCube_ARAT,
      spline_mod_ARAT)


# negative exponential model 
library(nlme)
set.seed(100)
negEx_mod_ARAT <- nlme(AffARATTotal ~ b_0i +
                           (b_1i)*(exp(b_2i * WeeksPostStroke)),
                         data = ACUTE,
                         fixed = b_0i + b_1i + b_2i ~ 1,
                         random = list(b_0i ~ 1, b_1i ~ 1, b_2i ~1),
                         groups = ~ SubIDName,
                         start = c(60, -60, -1),
                         na.action = na.omit)

summary(negEx_mod_ARAT)
-2 * as.numeric(logLik(negEx_mod_ARAT))
3467.8-3263.5
pchisq(q=3467.8-3263.5, df=3, lower.tail=FALSE)


# 3-parameter sigmoidal model
sigmoid_model <- function(x, a, b, x0) {
  a / (1 + exp(-b * (x - x0)))
}

sigmoid_mod_ARATmodel <- nlme(
  model = AffARATTotal ~ a / (1 + exp(-b * (WeeksPostStroke - x0))),  # Nonlinear model formula (sigmoid)
  data = ACUTE,                          # Your dataset
  fixed = a + b + x0 ~ 1,                    # Fixed effects for the parameters
  random = a + x0 ~ 1 | SubIDName,         # Random effects (allowing a, b, x0 to vary by subject)
  start = c(a = 60, b = 1, x0 = 5),         # Starting values for a, b, and x0
  na.action = na.omit,                       # Handling missing values by omitting rows with NAs
  control = nlmeControl(maxIter = 200, msMaxIter = 200)  # Control iterations for fitting
)

summary(sigmoid_mod_ARATmodel)
-2 * as.numeric(logLik(sigmoid_mod_ARATmodel))
3467.8-3263.5
pchisq(q=3467.8-3263.5, df=3, lower.tail=FALSE)



# Fitting Slopes for the Fugl Meyer ---------------------------------------------
# Intercept model
int_mod_FM<-lmer(UEFuglMeyer~
                # Fixed-effects
                1+
                # Random-effects
                (1|SubIDName), data=ACUTE, REML=FALSE,
              control=lmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=5e5)))
summary(int_mod_FM)

# linear model
lin_mod_FM<-lmer(UEFuglMeyer~
                # Fixed-effects
                1+
                  WeeksPostStroke+
                # Random-effects
                (1+WeeksPostStroke|SubIDName), data=ACUTE, REML=FALSE,
              control=lmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=5e5)))
summary(lin_mod_FM)

# quadratic model
sq_mod_FM<-lmer(UEFuglMeyer~
               # Fixed-effects
               1+WeeksPostStroke+I(WeeksPostStroke^2)+
               # Random-effects
               (1+WeeksPostStroke+I(WeeksPostStroke^2)|SubIDName), data=ACUTE, REML=FALSE,
             control=lmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=5e5)))
summary(sq_mod_FM)


# cubic model
cub_mod_FM<-lmer(UEFuglMeyer~
                # Fixed-effects
                1+WeeksPostStroke+I(WeeksPostStroke^2)+I(WeeksPostStroke^3)+
                # Random-effects
                (1+WeeksPostStroke+I(WeeksPostStroke^2)|SubIDName), data=ACUTE, REML=FALSE,
              control=lmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=5e5)))
summary(cub_mod_FM)

anova(int_mod_FM, lin_mod_FM, sq_mod_FM, cub_mod_FM)

# inverse square root
invSq_FM <-lmer(UEFuglMeyer~
               # Fixed-effects
               1+invTime.2+
               # Random-effects
               (1+invTime.2|SubIDName), data=ACUTE, REML=FALSE,
             control=lmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=5e5)))
summary(invSq_FM)

# inverse cube root
invCube_FM <-lmer(UEFuglMeyer~
                 # Fixed-effects
                 1+invTime.3+
                 # Random-effects
                 (1+invTime.3|SubIDName), data=ACUTE, REML=FALSE,
               control=lmerControl(optimizer="bobyqa",
                                   optCtrl=list(maxfun=5e5)))
summary(invCube_FM)
anova(int_mod_FM, lin_mod_FM, sq_mod_FM, cub_mod_FM, invSq_FM, invCube_FM)





# spline model
spline_mod_FM <-lmer(UEFuglMeyer~
                    # Fixed-effects
                    1+WeeksPostStroke+knot05+knot11+
                    # Random-effects
                    (1+WeeksPostStroke+knot05+knot11|SubIDName), data=ACUTE, REML=FALSE,
                  control=lmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=5e5)))
summary(spline_mod_FM)
anova(int_mod_FM, lin_mod_FM, sq_mod_FM, cub_mod_FM, invSq_FM, invCube_FM, spline_mod_FM)


# negative exponential model 
set.seed(100)
negEx_mod_FM <- nlme(UEFuglMeyer ~ b_0i +
                         (b_1i)*(exp(b_2i * WeeksPostStroke)),
                       data = ACUTE,
                       fixed = b_0i + b_1i + b_2i ~ 1,
                       random = list(b_0i ~ 1, b_1i ~ 1),
                       groups = ~ SubIDName,
                       start = c(70, -70, -1),
                       na.action = na.omit)

summary(negEx_mod_FM)
-2 * as.numeric(logLik(negEx_mod_FM))
2747.2-2564.951
pchisq(q=3467.8-3263.5, df=3, lower.tail=FALSE)







# Fitting Slopes for the Use Ratio ---------------------------------------------
# Intercept model
int_mod_UR<-lmer(use_ratio~
                # Fixed-effects
                1+
                # Random-effects
                (1|SubIDName), data=ACUTE, REML=FALSE,
              control=lmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=5e5)))
summary(int_mod_UR)

# linear model
lin_mod_UR<-lmer(use_ratio~
                # Fixed-effects
                1+WeeksPostStroke+
                # Random-effects
                (1+WeeksPostStroke|SubIDName), data=ACUTE, REML=FALSE,
              control=lmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=5e5)))
summary(lin_mod_UR)

# quadratic model
sq_mod_UR<-lmer(use_ratio~
               # Fixed-effects
               1+WeeksPostStroke+I(WeeksPostStroke^2)+
               # Random-effects
               (1+WeeksPostStroke+I(WeeksPostStroke^2)|SubIDName), data=ACUTE, REML=FALSE,
             control=lmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=5e5)))
summary(sq_mod_UR)


# cubic model
cub_mod_UR<-lmer(use_ratio~
                # Fixed-effects
                1+WeeksPostStroke+I(WeeksPostStroke^2)+I(WeeksPostStroke^3)+
                # Random-effects
                (1+WeeksPostStroke+I(WeeksPostStroke^2)|SubIDName), data=ACUTE, REML=FALSE,
              control=lmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=5e5)))
summary(cub_mod_UR)

anova(int_mod_UR, lin_mod_UR, sq_mod_UR, cub_mod_UR)

# inverse square root
invSq_UR <-lmer(use_ratio~
               # Fixed-effects
               1+invTime.2+
               # Random-effects
               (1+invTime.2|SubIDName), data=ACUTE, REML=FALSE,
             control=lmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=5e5)))
summary(invSq_UR)

# inverse cube root
invCube_UR <-lmer(use_ratio~
                 # Fixed-effects
                 1+invTime.3+
                 # Random-effects
                 (1+invTime.3|SubIDName), data=ACUTE, REML=FALSE,
               control=lmerControl(optimizer="bobyqa",
                                   optCtrl=list(maxfun=5e5)))
summary(invCube_UR)
anova(int_mod_UR, lin_mod_UR, sq_mod_UR, cub_mod_UR, invSq_UR, invCube_UR)





# spline model
spline_mod_UR <-lmer(use_ratio~
                    # Fixed-effects
                    1+WeeksPostStroke+knot05+knot11+
                    # Random-effects
                    (1+WeeksPostStroke+knot05+knot11|SubIDName), data=ACUTE, REML=FALSE,
                  control=lmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=5e5)))
summary(spline_mod_UR)
anova(int_mod_UR, lin_mod_UR, sq_mod_UR, cub_mod_UR, invSq_UR, invCube_UR, spline_mod_UR)




# negative exponential model 
set.seed(100)
negEx_mod_UR <- nlme(use_ratio ~ b_0i +
                         (b_1i)*(exp(b_2i * WeeksPostStroke)),
                       data = ACUTE,
                       fixed = b_0i + b_1i + b_2i ~ 1,
                       random = list(b_0i ~ 1),
                       groups = ~ SubIDName,
                       start = c(1, -1, -0.5),
                       na.action = na.omit,
                     control = lmeControl(maxIter = 500, msMaxIter = 500))

summary(negEx_mod_UR)
-2 * as.numeric(logLik(negEx_mod_UR))
-465.45--558.2836
pchisq(q=-465.45--558.2836, df=1, lower.tail=FALSE)




# 3-parameter sigmoidal model
sigmoid_model <- function(x, a, b, x0) {
  a / (1 + exp(-b * (x - x0)))
}

sigmoid_mod_URmodel <- nlme(
  model = use_ratio ~ a / (1 + exp(-b * (WeeksPostStroke - x0))),  # Nonlinear model formula (sigmoid)
  data = ACUTE,                          # Your dataset
  fixed = a + b + x0 ~ 1,                    # Fixed effects for the parameters
  random = a + x0 ~ 1 | SubIDName,         # Random effects (allowing a, b, x0 to vary by subject)
  start = c(a = 1.2, b = 1, x0 = 5),         # Starting values for a, b, and x0
  na.action = na.omit,                       # Handling missing values by omitting rows with NAs
  control = nlmeControl(maxIter = 500, msMaxIter = 500)  # Control iterations for fitting
)

summary(sigmoid_mod_URmodel)
-2 * as.numeric(logLik(sigmoid_mod_URmodel))
3467.8-3263.5
pchisq(q=3467.8-3263.5, df=3, lower.tail=FALSE)


# Extracting coefficients from each person from the best fitting models ---------
ARAT_COEFS <- coef(spline_mod_ARAT)$SubIDName %>% rownames_to_column(var="subID") %>%
  rename(atTime0 = `(Intercept)`,
         SlopeAt0 = WeeksPostStroke) %>%
  pivot_longer(cols=atTime0:knot11,
               names_to ="parameter",
               values_to="ARAT") 

FM_COEFS <- coef(spline_mod_FM)$SubIDName %>% rownames_to_column(var="subID") %>%
  rename(atTime0 = `(Intercept)`,
         SlopeAt0 = WeeksPostStroke) %>%
  pivot_longer(cols=atTime0:knot11,
               names_to ="parameter",
               values_to="FM") 

UR_COEFS <- coef(spline_mod_UR)$SubIDName %>% rownames_to_column(var="subID") %>%
  rename(atTime0 = `(Intercept)`,
         SlopeAt0 = WeeksPostStroke) %>%
  pivot_longer(cols=atTime0:knot11,
               names_to ="parameter",
               values_to="UR") 

COEFS <- merge(x=ARAT_COEFS,
               y=FM_COEFS,
               by=c("subID", "parameter"), all=TRUE)

COEFS <- merge(x=COEFS,
               y=UR_COEFS,
               by=c("subID", "parameter"), all=TRUE)

head(COEFS)




# Mutate the coefficients to wide format
# Calculate series of slopes and intercepts FROM SPLINE KNOTS ------------------
# Although knots are directly available from the regression, transforming these
# knots into a series of slopes and intercepts appears to be more interpretable.
COEFS_WIDE <- COEFS %>% pivot_wider(values_from=ARAT:UR,
                               names_from = parameter,
                               names_sep = "_") %>%
  mutate(ARAT_atTime05=ARAT_atTime0+ARAT_SlopeAt0*5,
         ARAT_SlopeAt05=ARAT_SlopeAt0+ARAT_knot05,
         ARAT_atTime11=ARAT_atTime0+ARAT_SlopeAt0*11+ARAT_knot05*6,
         ARAT_SlopeAt11=ARAT_SlopeAt0+ARAT_knot05+ARAT_knot11,
         #
         FM_atTime05=FM_atTime0+FM_SlopeAt0*5,
         FM_SlopeAt05=FM_SlopeAt0+FM_knot05,
         FM_atTime11=FM_atTime0+FM_SlopeAt0*11+FM_knot05*6,
         FM_SlopeAt11=FM_SlopeAt0+FM_knot05+FM_knot11,
         #
         UR_atTime05=UR_atTime0+UR_SlopeAt0*5,
         UR_SlopeAt05=UR_SlopeAt0+UR_knot05,
         UR_atTime11=UR_atTime0+UR_SlopeAt0*11+UR_knot05*6,
         UR_SlopeAt11=UR_SlopeAt0+UR_knot05+UR_knot11
  )

head(COEFS)
head(COEFS_WIDE)

# correlations between coefficients --------------------------------------------
COEFS_LONG <- COEFS_WIDE %>%
  select(subID, 
         ARAT_atTime0, FM_atTime0, UR_atTime0,
         ARAT_SlopeAt0, FM_SlopeAt0, UR_SlopeAt0,
         ARAT_SlopeAt05, FM_SlopeAt05, UR_SlopeAt05,
         ARAT_SlopeAt11, FM_SlopeAt11, UR_SlopeAt11) %>%
  pivot_longer(cols=ARAT_atTime0:UR_SlopeAt11,
               names_to = c("Variable", "Time"),
               names_sep = "_",
               values_to = "Coefficient") %>%
  pivot_wider(values_from = "Coefficient", 
              names_from = "Variable")

head(COEFS_LONG)

time_list <- c("atTime0", "SlopeAt0", "SlopeAt05", "SlopeAt11")

for(t in time_list) {
  print(t)
  
  print("ARAT -- FM UE")
  print(cor.test(x=COEFS_LONG[COEFS_LONG$Time==eval(t),]$ARAT,
           y=COEFS_LONG[COEFS_LONG$Time==eval(t),]$FM),
        method="pearson",
        conf.level = 0.95)
  
  print("Use Ratio -- ARAT")
  print(cor.test(x=COEFS_LONG[COEFS_LONG$Time==eval(t),]$UR,
                 y=COEFS_LONG[COEFS_LONG$Time==eval(t),]$ARAT),
        method="pearson",
        conf.level = 0.95)
  
  print("Use Ratio -- FM")
  print(cor.test(x=COEFS_LONG[COEFS_LONG$Time==eval(t),]$FM,
                 y=COEFS_LONG[COEFS_LONG$Time==eval(t),]$UR),
        method="pearson",
        conf.level = 0.95)
  
}


# Figure 4: Spaghetti plots of trajectories ----
# PMC8442937_049
# PMC8442937_017
colnames(ACUTE)

COEFS_WIDE %>% filter(subID=="PMC8442937_017")

ggplot(data=ACUTE %>%
         filter(SubIDName=="PMC8442937_017"),
       aes(x=WeeksPostStroke, y=AffARATTotal))+
  geom_line(aes(group=SubIDName), col="black", lty=1,
            position = position_dodge(width=0.4))+
  geom_point(aes(group=SubIDName), col="black", shape=21,
             position = position_dodge(width=0.4))+
  geom_segment(aes(x = 0, 
                   y = -3.54+1.21*0+0.660*0-1.44*0, 
                   xend = 5, 
                   yend = -3.54+1.21*5+0.660*0-1.44*0),
               col=cbPalette[3],
               arrow = arrow(length = unit(0.2,"cm")), lwd=1.5)+
  geom_segment(aes(x = 5, 
                   y = -3.54+1.21*5+0.660*0-1.44*0, 
                   xend = 11, 
                   yend = -3.54+1.21*11+0.660*6-1.44*0),
               col=cbPalette[4],
               arrow = arrow(length = unit(0.2,"cm")), lwd=1.5)+
  geom_segment(aes(x = 11, 
                   y = -3.54+1.21*11+0.660*6-1.44*0, 
                   xend = 25, 
                   yend = -3.54+1.21*25+0.660*20-1.44*14),
               col=cbPalette[2],
               arrow = arrow(length = unit(0.2,"cm")), lwd=1.5)+
  scale_x_continuous(name = "Weeks Post-Stroke", limits=c(0,36)) +
  scale_y_continuous(name = "ARAT (Affected Arm)", limits=c(-5,60),
                     breaks = c(seq(0,60,10))) +
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/ARAT_splne_singleSubject_A.jpeg",
  plot = last_plot(),
  width = 3,
  height = 2.5,
  units = "in",
  dpi = 300
)



COEFS_WIDE %>% filter(subID=="PMC8442937_049")

ggplot(data=ACUTE %>%
         filter(SubIDName=="PMC8442937_049"),
       aes(x=WeeksPostStroke, y=AffARATTotal))+
  geom_line(aes(group=SubIDName), col="black", lty=1,
            position = position_dodge(width=0.4))+
  geom_point(aes(group=SubIDName), col="black", shape=21,
             position = position_dodge(width=0.4))+
  geom_segment(aes(x = 0, 
                   y = -26.0+10.0*0-5.42*0-4.66*0, 
                   xend = 5, 
                   yend = -26.0+10.0*5-5.42*0-4.66*0),
               col=cbPalette[3],
               arrow = arrow(length = unit(0.2,"cm")), lwd=1.5)+
  geom_segment(aes(x = 5, 
                   y = -26.0+10.0*5-5.42*0-4.66*0, 
                   xend = 11, 
                   yend = -26.0+10.0*11-5.42*6-4.66*0),
               col=cbPalette[4],
               arrow = arrow(length = unit(0.2,"cm")), lwd=1.5)+
  geom_segment(aes(x = 11, 
                   y = -26.0+10.0*11-5.42*6-4.66*0, 
                   xend = 22, 
                   yend = -26.0+10.0*22-5.42*17-4.66*11),
               col=cbPalette[2],
               arrow = arrow(length = unit(0.2,"cm")), lwd=1.5)+
  scale_x_continuous(name = "Weeks Post-Stroke", limits=c(0,36)) +
  scale_y_continuous(name = "ARAT (Affected Arm)", limits=c(-30,60),
                     breaks = c(seq(0,60,10))) +
  #facet_wrap(~enrollTime, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/ARAT_splne_singleSubject_B.jpeg",
  plot = last_plot(),
  width = 3,
  height = 2.5,
  units = "in",
  dpi = 300
)

# Figure 4A : ARAT Splines ----
head(COEFS_WIDE)
head(UNIQUE)
COEFS_WIDE <- merge(x=COEFS_WIDE, 
                    y=UNIQUE %>% select(SubIDName, StrokeType),
                    by.x="subID",
                    by.y="SubIDName")

ggplot(data=COEFS_WIDE)+
  geom_segment(aes(group=subID,
                   x = 0, 
                   y = ARAT_atTime0, 
                   xend = 5, 
                   yend = ARAT_atTime0+ARAT_SlopeAt0*5),
               col=cbPalette[3],
               arrow = arrow(length = unit(0.2,"cm")), lwd=0.5, alpha=0.5)+
  geom_segment(aes(group=subID,
                   x = 5, 
                   y = ARAT_atTime0+ARAT_SlopeAt0*5, 
                   xend = 11, 
                   yend = ARAT_atTime0+ARAT_SlopeAt0*11+ARAT_knot05*6),
               col=cbPalette[4],
               arrow = arrow(length = unit(0.2,"cm")), lwd=0.5, alpha=0.5)+
  geom_segment(aes(group=subID,
                   x = 11, 
                   y = ARAT_atTime0+ARAT_SlopeAt0*11+ARAT_knot05*6, 
                   xend = 23, 
                   yend = ARAT_atTime0+ARAT_SlopeAt0*23+ARAT_knot05*18+ARAT_knot11*12),
               col=cbPalette[2],
               arrow = arrow(length = unit(0.2,"cm")), lwd=0.5, alpha=0.5)+
  scale_x_continuous(name = "Weeks Post-Stroke", limits=c(0,24)) +
  scale_y_continuous(name = "ARAT (Affected Arm)", limits=c(-30,60),
                     breaks = c(seq(0,60,10))) +
  #facet_wrap(~StrokeType, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/ARAT_splines.jpeg",
  plot = last_plot(),
  width = 4,
  height = 2.5,
  units = "in",
  dpi = 300
)


# Figure 4B : FM Splines ----
head(COEFS_WIDE)
ggplot(data=COEFS_WIDE)+
  geom_segment(aes(group=subID,
                   x = 0, 
                   y = FM_atTime0, 
                   xend = 5, 
                   yend = FM_atTime0+FM_SlopeAt0*5),
               col=cbPalette[3],
               arrow = arrow(length = unit(0.2,"cm")), lwd=0.5, alpha=0.5)+
  geom_segment(aes(group=subID,
                   x = 5, 
                   y = FM_atTime0+FM_SlopeAt0*5, 
                   xend = 11, 
                   yend = FM_atTime0+FM_SlopeAt0*11+FM_knot05*6),
               col=cbPalette[4],
               arrow = arrow(length = unit(0.2,"cm")), lwd=0.5, alpha=0.5)+
  geom_segment(aes(group=subID,
                   x = 11, 
                   y = FM_atTime0+FM_SlopeAt0*11+FM_knot05*6, 
                   xend = 23, 
                   yend = FM_atTime0+FM_SlopeAt0*23+FM_knot05*18+FM_knot11*12),
               col=cbPalette[2],
               arrow = arrow(length = unit(0.2,"cm")), lwd=0.5, alpha=0.5)+
  scale_x_continuous(name = "Weeks Post-Stroke", limits=c(0,24)) +
  scale_y_continuous(name = "FM UE", limits=c(-5,70),
                     breaks = c(seq(0,70,10))) +
  #facet_wrap(~StrokeType, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/FM_splines.jpeg",
  plot = last_plot(),
  width = 4,
  height = 2.5,
  units = "in",
  dpi = 300
)




# Figure 4C : UR Splines ----
head(COEFS_WIDE)
ggplot(data=COEFS_WIDE)+
  geom_segment(aes(group=subID,
                   x = 0, 
                   y = UR_atTime0, 
                   xend = 5, 
                   yend = UR_atTime0+UR_SlopeAt0*5),
               col=cbPalette[3],
               arrow = arrow(length = unit(0.2,"cm")), lwd=0.5, alpha=0.5)+
  geom_segment(aes(group=subID,
                   x = 5, 
                   y = UR_atTime0+UR_SlopeAt0*5, 
                   xend = 11, 
                   yend = UR_atTime0+UR_SlopeAt0*11+UR_knot05*6),
               col=cbPalette[4],
               arrow = arrow(length = unit(0.2,"cm")), lwd=0.5, alpha=0.5)+
  geom_segment(aes(group=subID,
                   x = 11, 
                   y = UR_atTime0+UR_SlopeAt0*11+UR_knot05*6, 
                   xend = 23, 
                   yend = UR_atTime0+UR_SlopeAt0*23+UR_knot05*18+UR_knot11*12),
               col=cbPalette[2],
               arrow = arrow(length = unit(0.2,"cm")), lwd=0.5, alpha=0.5)+
  scale_x_continuous(name = "Weeks Post-Stroke", limits=c(0,24)) +
  scale_y_continuous(name = "Use Ratio", limits=c(-0.2,1.4),
                     breaks = c(seq(0,1.4,0.2))) +
  #facet_wrap(~StrokeType, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(size=12, color="black"),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")

ggsave(
  filename="./outputs/UR_splines.jpeg",
  plot = last_plot(),
  width = 4,
  height = 2.5,
  units = "in",
  dpi = 300
)





# Figure 5: Density and covariance of different parameters ----
colnames(COEFS_WIDE)

# format control
panel_settings <- 'axis.text=element_text(size=0, color="black"),
        legend.text=element_text(size=8, color="black"),
        legend.title=element_text(size=8, face="bold"),
        axis.title=element_text(size=8, face="bold"),
        plot.title=element_text(size=8, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=8, face="bold"),
        legend.position = "none"'

# make a simple function to scale x-axis values
scaleFUN <- function(x) sprintf("%.2f", x) # 1f = 1 decimal, 2f=two decimals, etc
# This is only really needed on use ratio slope plots



ARAT_int <- ggplot(data=COEFS_WIDE,
       aes(x=ARAT_atTime0))+
  geom_histogram(fill=cbPalette[3], col="black", lty=1,
            position = position_dodge(width=0.4))+
  scale_x_continuous(name = "Est. ARAT") +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

FM_int <- ggplot(data=COEFS_WIDE,
       aes(x=FM_atTime0))+
  geom_histogram(fill=cbPalette[3], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "Est. FM UE") +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

UR_int <- ggplot(data=COEFS_WIDE,
       aes(x=UR_atTime0))+
  geom_histogram(fill=cbPalette[3], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "Est. Use Ratio") +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

UR_FM <- ggplot(data=COEFS_WIDE,
       aes(x=UR_atTime0,
           y=FM_atTime0))+
  geom_point(fill=cbPalette[3], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "Est. Use Ratio") +
  scale_y_continuous(name = "Est. FM UE") +
  theme_bw()+
  theme(eval(panel_settings))

UR_ARAT <- ggplot(data=COEFS_WIDE,
                aes(x=UR_atTime0,
                    y=ARAT_atTime0))+
  geom_point(fill=cbPalette[3], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "Est. Use Ratio") +
  scale_y_continuous(name = "Est. ARAT") +
  theme_bw()+
  theme(eval(panel_settings))



FM_ARAT <- ggplot(data=COEFS_WIDE,
                  aes(x=FM_atTime0,
                      y=ARAT_atTime0))+
  geom_point(fill=cbPalette[3], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "Est. FM UE") +
  scale_y_continuous(name = "Est. ARAT") +
  theme_bw()+
  theme(eval(panel_settings))


(ARAT_int + plot_spacer() + plot_spacer())/
  (FM_ARAT + FM_int + plot_spacer())/
  (UR_ARAT + UR_FM + UR_int)

ggsave(
  filename="./outputs/w0_intercepts.jpeg",
  plot = last_plot(),
  width = 6,
  height = 5,
  units = "in",
  dpi = 150
)


# ARAT Slopes over time ----
ARAT00 <- ggplot(data=COEFS_WIDE,
       aes(x=ARAT_SlopeAt0))+
  geom_histogram(fill=cbPalette[3], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "ARAT Slopes", limits=c(-2, 15)) +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

ARAT05 <- ggplot(data=COEFS_WIDE,
       aes(x=ARAT_SlopeAt05))+
  geom_histogram(fill=cbPalette[4], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "ARAT Slopes", limits=c(-2, 15)) +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

ARAT11 <- ggplot(data=COEFS_WIDE,
       aes(x=ARAT_SlopeAt11))+
  geom_histogram(fill=cbPalette[2], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "ARAT Slopes", limits=c(-2, 15)) +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

ARAT00/ARAT05/ARAT11

ggsave(
  filename="./outputs/arat_slopes.jpeg",
  plot = last_plot(),
  width = 2.5,
  height = 6,
  units = "in",
  dpi = 150
)



# FM Slopes over time ----
FM00 <- ggplot(data=COEFS_WIDE,
               aes(x=FM_SlopeAt0))+
  geom_histogram(fill=cbPalette[3], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "FMUE Slopes", limits=c(-1, 8)) +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

FM05 <- ggplot(data=COEFS_WIDE,
               aes(x=FM_SlopeAt05))+
  geom_histogram(fill=cbPalette[4], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "FMUE Slopes", limits=c(-1, 8)) +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))


FM11 <- ggplot(data=COEFS_WIDE,
               aes(x=FM_SlopeAt11))+
  geom_histogram(fill=cbPalette[2], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "FMUE Slopes", limits=c(-1, 8)) +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

FM00/FM05/FM11

ggsave(
  filename="./outputs/fmue_slopes.jpeg",
  plot = last_plot(),
  width = 2.5,
  height = 6,
  units = "in",
  dpi = 150
)






# Use Ratio Slopes over time ----
UR00 <- ggplot(data=COEFS_WIDE,
               aes(x=UR_SlopeAt0))+
  geom_histogram(fill=cbPalette[3], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "Use Ratio Slopes", limits=c(-0.05,0.15)) +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

UR05 <- ggplot(data=COEFS_WIDE,
               aes(x=UR_SlopeAt05))+
  geom_histogram(fill=cbPalette[4], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "Use Ratio Slopes", limits=c(-0.05,0.15)) +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

UR11 <- ggplot(data=COEFS_WIDE,
               aes(x=UR_SlopeAt11))+
  geom_histogram(fill=cbPalette[2], col="black", lty=1,
                 position = position_dodge(width=0.4))+
  scale_x_continuous(name = "Use Ratio Slopes", limits=c(-0.05,0.15)) +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  theme(eval(panel_settings))

UR00/UR05/UR11

ggsave(
  filename="./outputs/useRatio_slopes.jpeg",
  plot = last_plot(),
  width = 2.5,
  height = 6,
  units = "in",
  dpi = 150
)


# Week 0 Slopes ----------------------------------------------------------------
colnames(COEFS_WIDE)
UR_FM_00 <- ggplot(data=COEFS_WIDE,
                aes(x=UR_SlopeAt0,
                    y=FM_SlopeAt0))+
  geom_point(fill=cbPalette[3], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "Use Ratio Slope") +
  scale_y_continuous(name = "FM UE Slope") +
  theme_bw()+
  theme(eval(panel_settings))

UR_ARAT_00 <- ggplot(data=COEFS_WIDE,
                     aes(x=UR_SlopeAt0,
                         y=ARAT_SlopeAt0))+
  geom_point(fill=cbPalette[3], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "Use Ratio Slope") +
  scale_y_continuous(name = "ARAT Slope") +
  theme_bw()+
  theme(eval(panel_settings))



FM_ARAT_00 <- ggplot(data=COEFS_WIDE,
                     aes(x=FM_SlopeAt0,
                         y=ARAT_SlopeAt0))+
  geom_point(fill=cbPalette[3], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "FM UE Slope") +
  scale_y_continuous(name = "ARAT Slope") +
  theme_bw()+
  theme(eval(panel_settings))


(ARAT00 + plot_spacer() + plot_spacer())/
  (FM_ARAT_00 + FM00 + plot_spacer())/
  (UR_ARAT_00 + UR_FM_00 + UR00)

ggsave(
  filename="./outputs/w0_slopes.jpeg",
  plot = last_plot(),
  width = 6,
  height = 5,
  units = "in",
  dpi = 150
)




# Week 5 Slopes ----------------------------------------------------------------
colnames(COEFS_WIDE)
UR_FM_05 <- ggplot(data=COEFS_WIDE,
                   aes(x=UR_SlopeAt05,
                       y=FM_SlopeAt05))+
  geom_point(fill=cbPalette[4], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "Use Ratio Slope", labels=scaleFUN) +
  scale_y_continuous(name = "FM UE Slope") +
  theme_bw()+
  theme(eval(panel_settings))

UR_ARAT_05 <- ggplot(data=COEFS_WIDE,
                     aes(x=UR_SlopeAt05,
                         y=ARAT_SlopeAt05))+
  geom_point(fill=cbPalette[4], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "Use Ratio Slope", labels=scaleFUN) +
  scale_y_continuous(name = "ARAT Slope") +
  theme_bw()+
  theme(eval(panel_settings))



FM_ARAT_05 <- ggplot(data=COEFS_WIDE,
                     aes(x=FM_SlopeAt05,
                         y=ARAT_SlopeAt05))+
  geom_point(fill=cbPalette[4], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "FM UE Slope") +
  scale_y_continuous(name = "ARAT Slope") +
  theme_bw()+
  theme(eval(panel_settings))


(ARAT05 + plot_spacer() + plot_spacer())/
  (FM_ARAT_05 + FM05 + plot_spacer())/
  (UR_ARAT_05 + UR_FM_05 + UR05)

ggsave(
  filename="./outputs/w05_slopes.jpeg",
  plot = last_plot(),
  width = 6,
  height = 5,
  units = "in",
  dpi = 150
)




# Week 11 Slopes ----------------------------------------------------------------
colnames(COEFS_WIDE)
UR_FM_11 <- ggplot(data=COEFS_WIDE,
                   aes(x=UR_SlopeAt11,
                       y=FM_SlopeAt11))+
  geom_point(fill=cbPalette[2], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "Use Ratio Slope", labels=scaleFUN) +
  scale_y_continuous(name = "FM UE Slope") +
  theme_bw()+
  theme(eval(panel_settings))

UR_ARAT_11 <- ggplot(data=COEFS_WIDE,
                     aes(x=UR_SlopeAt11,
                         y=ARAT_SlopeAt11))+
  geom_point(fill=cbPalette[2], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "Use Ratio Slope", labels=scaleFUN) +
  scale_y_continuous(name = "ARAT Slope") +
  theme_bw()+
  theme(eval(panel_settings))



FM_ARAT_11 <- ggplot(data=COEFS_WIDE,
                     aes(x=FM_SlopeAt11,
                         y=ARAT_SlopeAt11))+
  geom_point(fill=cbPalette[2], shape=21, alpha=0.8)+
  stat_poly_line() +
  stat_poly_eq() +
  scale_x_continuous(name = "FM UE Slope") +
  scale_y_continuous(name = "ARAT Slope") +
  theme_bw()+
  theme(eval(panel_settings))


(ARAT11 + plot_spacer() + plot_spacer())/
  (FM_ARAT_11 + FM11 + plot_spacer())/
  (UR_ARAT_11 + UR_FM_11 + UR11)

ggsave(
  filename="./outputs/w11_slopes.jpeg",
  plot = last_plot(),
  width = 6,
  height = 5,
  units = "in",
  dpi = 150
)






# Correlations among parameters from the same model ----------
colnames(COEFS_WIDE)
ggpairs(COEFS_WIDE %>% select(-subID),
        columns = c("ARAT_atTime0", "ARAT_SlopeAt0",
                    "ARAT_SlopeAt05", "ARAT_SlopeAt11"),
        mapping = ggplot2::aes(alpha=0.5),
        lower = list(continuous = wrap("points", shape=21)),
        upper = list(continuous = wrap("cor", size = 5)),
        diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag"))+
  theme_bw()+scale_color_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)

ggsave(
  filename="./outputs/correlations_ARAT.jpeg",
  plot = last_plot(),
  width = 6,
  height = 5,
  units = "in",
  dpi = 150
)


ggpairs(COEFS_WIDE %>% select(-subID),
        columns = c("FM_atTime0", "FM_SlopeAt0",
                    "FM_SlopeAt05", "FM_SlopeAt11"),
        mapping = ggplot2::aes(alpha=0.5),
        lower = list(continuous = wrap("points", shape=21)),
        upper = list(continuous = wrap("cor", size = 5)),
        diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag"))+
  theme_bw()+scale_color_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)

ggsave(
  filename="./outputs/correlations_FM.jpeg",
  plot = last_plot(),
  width = 6,
  height = 5,
  units = "in",
  dpi = 150
)




ggpairs(COEFS_WIDE %>% select(-subID),
        columns = c("UR_atTime0", "UR_SlopeAt0",
                    "UR_SlopeAt05", "UR_SlopeAt11"),
        mapping = ggplot2::aes(alpha=0.5),
        lower = list(continuous = wrap("points", shape=21)),
        upper = list(continuous = wrap("cor", size = 5)),
        diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag"))+
  theme_bw()+scale_color_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)

ggsave(
  filename="./outputs/correlations_UR.jpeg",
  plot = last_plot(),
  width = 6,
  height = 5,
  units = "in",
  dpi = 150
)





# Part 2: Demographic statistics for the Longitudinal cohort --------------------------------
# Confirm correct filtering of acute longitudinal cohort
summary(LONG$enrollCat)
summary(ACUTE$enrollCat)
# extract subject IDs that we want
sublist <- c(unique(ACUTE$SubIDName))

# Use the same enrollment data as before, but now subset on the longitudinal IDs
UNIQUE_ACUTE<-UNIQUE %>% filter(SubIDName %in% sublist)

colnames(UNIQUE_ACUTE)

UNIQUE_ACUTE %>% summary()
summary(as.factor(UNIQUE_ACUTE$NumberofStrokes))





