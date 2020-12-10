library(tidyverse)
library(lmerTest)
library(emmeans)

gait_nl_raw = read.csv('/Volumes/klouie1/Data_Analysis/No_level_gait_metrics.csv')
gait_sl_raw = read.csv('/Volumes/klouie1/Data_Analysis/Single_level_gait_metrics.csv')

gait_nl = gait_nl_raw %>% filter(Step_Time < 1.5, Step_Length < 1.5, Stride_Time < 1.5, Stride_Length < 1.5) %>% mutate(Level_Type = rep('None',n()))
gait_sl = gait_sl_raw %>% filter(Step_Time < 1.5, Step_Length < 1.5, Stride_Time < 1.5, Stride_Length < 1.5) %>% mutate(Level_Type = rep('Single',n()))
gait = gait_nl %>% bind_rows(gait_sl)

gait_avg = gait %>% 
  group_by(Side,Level_Type) %>%
  summarise_all(list(~mean(.,na.rm = T),~sd(.,na.rm = T)))

ggplot(data = gait_avg,aes(x = Side, y = Step_Time_mean, fill = Level_Type)) +
  geom_col(position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin = Step_Time_mean - Step_Time_sd, ymax = Step_Time_mean + Step_Time_sd), width = 0.2, position = position_dodge(width = 1)) + 
  ylab('Step Time (s)') +
  theme_classic()

ggplot(data = gait_avg,aes(x = Side, y = Step_Length_mean, fill = Level_Type)) +
  geom_col(position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin = Step_Length_mean - Step_Length_sd, ymax = Step_Length_mean + Step_Length_sd), width = 0.2, position = position_dodge(width = 1)) + 
  ylab('Step Length (m)') +
  theme_classic()

ggplot(data = gait_avg,aes(x = Side, y = Step_Width_mean, fill = Level_Type)) +
  geom_col(position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin = Step_Width_mean - Step_Width_sd, ymax = Step_Width_mean + Step_Width_sd), width = 0.2, position = position_dodge(width = 1)) + 
  ylab('Step Width (m)') +
  theme_classic()
  
ggplot(data = gait_avg,aes(x = Side, y = Stride_Time_mean, fill = Level_Type)) +
  geom_col(position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin = Stride_Time_mean - Stride_Time_sd, ymax = Stride_Time_mean + Stride_Time_sd), width = 0.2, position = position_dodge(width = 1)) + 
  ylab('Stride Time (s)') +
  theme_classic()

ggplot(data = gait_avg,aes(x = Side, y = Stride_Length_mean, fill = Level_Type)) +
  geom_col(position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin = Stride_Length_mean - Stride_Length_sd, ymax = Stride_Length_mean + Stride_Length_sd), width = 0.2, position = position_dodge(width = 1)) + 
  ylab('Stride Length (m)') +
  theme_classic()

