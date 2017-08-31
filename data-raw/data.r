log <- file("r.log")
sink(log)
sink(log, type="message")
library(knitr)

render_listings()
knit("data.rnw", "data.tex")
sink()
sink(type="message")