library(mdtrt)

# data source
library(tswbench)

# data.table
library(data.table)

# warm-up data source
codes <- tswbench::default_srt_codes(stock = TRUE, fund = FALSE, index = FALSE)
dt <- tswbench::sina_realtime_quote(codes)

# make reused data matrix
n <- length(codes)
m <- matrix(0.0, nrow = n, ncol = 8)
colnames(m) <- c("Open", "High", "Low", "Close", "Tnvr", "Vol", "VWAP", "Time")

# make sliding window clock

# sliding window size
window_size_second <- 15
# provide enough buffer to store window data, 16 is more than enough since sina provides 3s snapshots
buffer_size <- 16
clk <- new(mdtrt:::sliding_window_clock_market, n, 15, 16)

while (TRUE) {

  dt <- tswbench::sina_realtime_quote(codes)
  # ensure row order is consistent, not necessary when using sina_realtime_quote()
  data.table::setkeyv(dt, "sina_code")

  # update_tvol() since dt$Vol is total volume traded, if Vol is snapshot volume, use update() instead
  clk$update_tvol(dt$Time, dt$Price, dt$Vol, m)

  # notice that m is updated by-ref
  # do stuff with m
  hist(m[, "Tnvr"])

  # stateful online analysis recommended, stay tune for updates in a few days
  # ema <- ewcov$update(m[, "VWAP"])

  #sleep(2.0)
}
