library(fpp3)
library(here)
d <- read.csv(here("obj/chunks/chunk_500.csv"))
d <- d[,-1]
a <- d[1,]
a <- a |> unlist(x = _)
attributes(a) <- NULL

series <- ts(a/10000,start = c(1981,7), frequency = 12)
plot(series)
series <- as_tsibble(series)
autoplot(series) + theme_minimal()
series <- fill_gaps(series)
saveRDS(series,here("obj/series.Rds"))
dcmp <- series %>%
  model(stl = STL(value))
components(dcmp)
components(dcmp) %>% autoplot()
ggsave("series.pdf")
