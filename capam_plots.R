#-----------------------------------------------------------------------------------------------------
#Load the Data
the_data <- process_bigeye_data(nspp = 3)

the_data[[2]] %>% group_by(lat, lon, year)

#Plot the data for capam 
world_map <- map_data("world")  

####Map of average counts data
counts <- the_data[[1]] %>% group_by(lat, lon) %>% summarize(min_cpue = min(cpue),
  avg_cpue = mean(cpue), max_cpue = max(cpue), diff_cpue = max_cpue - min_cpue,
  sd_cpue = sd(cpue), cv_cpue = sd_cpue / avg_cpue) %>% as.data.frame

ggplot(counts, aes(x = lon, y = lat)) + geom_tile(aes(fill = avg_cpue)) +
  geom_map(data = world_map, map = world_map, 
    aes(x = long, y = lat, map_id = region)) + 
  scale_x_continuous(limits = c(-148, -72)) + 
  scale_y_continuous(limits = c(-43, 48)) + 
  scale_fill_gradient(low = 'white', high = 'red') + theme_bw() + 
  ggsave("bigeye_counts.png", width = 7, height = 7)

####Map of average counts data
ggplot(counts, aes(x = lon, y = lat)) + geom_tile(aes(fill = cv_cpue)) +
  geom_map(data = world_map, map = world_map, 
    aes(x = long, y = lat, map_id = region)) + 
  scale_x_continuous(limits = c(-148, -72)) + 
  scale_y_continuous(limits = c(-43, 48)) + 
  scale_fill_gradient(low = 'white', high = 'blue') + theme_bw() +
  ggsave("bigeye_cv.png", width = 7, height = 7)

####Map of past 20 years
cpue_hook <- the_data[[1]]
cpue_hook$cpue <- cpue_hook$cpue / cpue_hook$hooks
cpue_hook[which(is.na(cpue_hook$cpue)), 'cpue'] <- 0

cpue_hook %>% filter(year >= 1990) %>% ggplot(aes(x = lon, y = lat)) +
  geom_tile(aes(fill = cpue)) +
  geom_map(data = world_map, map = world_map, 
    aes(x = long, y = lat, map_id = region)) + 
  scale_x_continuous(limits = c(-148, -72)) + 
  scale_y_continuous(limits = c(-43, 48)) + 
  scale_fill_gradient(low = 'white', high = 'red') + theme_bw() + 
  facet_wrap(~ year) + 
  ggsave('bigeye_cpue.png', width = 9.5, height = 8.8)

#--------------------------------------------------------------
####Map of japanese catch compositions
jpn_cc <- the_data[[2]]

jpn_cc %>% group_by(lat, lon, spp) %>% summarize(avg_cpue = mean(cpue)) %>% 
  ggplot(aes(x = lon, y = lat)) + geom_tile(aes(fill = avg_cpue)) + 
  facet_wrap(~ spp) + scale_fill_gradient(low = 'white', high = 'red') +
  geom_map(data = world_map, map = world_map, 
    aes(x = long, y = lat, map_id = region)) + 
  scale_x_continuous(limits = c(-148, -72)) + 
  scale_y_continuous(limits = c(-43, 48)) + theme_bw() + 
  ggsave("bigeye_size_counts.png", width = 12, height = 6)

#Look at size classes
jpn_cc %>% group_by(spp) %>% distinct(bin_range) %>% filter(spp == 3)

