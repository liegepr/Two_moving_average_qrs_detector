require("plotly")
source("qrs_detector_2ma_mod.r")

ad100 <- fread("100.csv", sep = ",")
                      
# ==================== static diagnostic plots ==========================================================
ad100_qrs <- ma_detector(signal = ad100[, MLII], srate = 360L, lowcut_f1 = 8L, highcut_f2 = 21L, filter_order = 3L, 
     qrs_win1 = 35L, beat_win2 = 220L, srate_ref = 360L, offset = 0.08, offset_win3 = 10L,
    slackness_red = TRUE, slackness_win1 = 0.2, slackness_win2 = 0.14, refractory_period = 0.3)

ad100_qrs0 <- ad100_qrs[648600:650000]
ad100_qrs0[!is.na(block.clean), block.clean := 1]
ad100_qrs0[is.na(block.clean), block.clean := 0]

# plot the squared signal along with windows, blocks and estimated locations
plot(signal_squared ~ idx, data = ad100_qrs0, type = "l", xlab = "sample (360 Hz)")
lines(mwa_qrs ~ idx, data = ad100_qrs0, type = "l", col = 2) 
lines(mwa_beat ~ idx, data = ad100_qrs0, type = "l", col = 3) 
lines(I(0.2 * block.clean) ~ idx, data = ad100_qrs0, type = "l", col = 5) 
points(signal_squared ~ loc, data = ad100_qrs0, col = "darkred")
points(signal_squared ~ loc_sr, data = ad100_qrs0, col = "#036830", pch = 3)
legend(x = "topleft",          # Position
       title = "Slackness reduction",
       legend = c("No", "Yes", "W1", "W2", "Block"),  # Legend texts
       pch = c(1, 3, 0, 0, 0),           # Line types
       col = c("darkred", "#036830", 3, 2, 5),           # Line colors
       lty = c(0, 0, 1, 1, 1))                 # Line width

# same as above, using time as x-axis
plot(signal_squared ~ time_stamp, data = ad100_qrs0, type = "l", xlab = "time (s)")
lines(mwa_qrs ~ time_stamp, data = ad100_qrs0, type = "l", col = 2) 
lines(mwa_beat ~ time_stamp, data = ad100_qrs0, type = "l", col = 3) 
lines(I(0.2 * block.clean) ~ time_stamp, data = ad100_qrs0, type = "l", col = 5) 
points(signal_squared ~ I(loc / 360), data = ad100_qrs0, col = "darkred")
points(signal_squared ~ I(loc / 360), data = ad100_qrs0, col = "#036830", pch = 3)
legend(x = "topleft",          # Position
       title = "Slackness reduction",
       legend = c("No", "Yes", "W1", "W2", "Block"),  # Legend texts
       pch = c(1, 3, 0, 0, 0),           # Line types
       col = c("darkred", "#036830", 3, 2, 5),           # Line colors
       lty = c(0, 0, 1, 1, 1))                 # Line width

# plot the unfiltered signal and check location accuracy
ad100_qrs1 <- ad100_qrs[649500:650000]
plot(signal ~ idx, data = ad100_qrs1, type = "l", xlab = "sample (360 Hz)")
points(signal ~ loc, data = ad100_qrs1, col = "darkred")
points(signal ~ loc_sr, data = ad100_qrs1, col = "#036830", pch = 3)
legend(x = "topleft",          # Position
       title = "Slackness reduction",
       legend = c("No", "Yes"),  # Legend texts
       pch = c(1, 3),           # Line types
       col = c("darkred", "#036830"),           # Line colors
       lty = 0)                 # Line width

# -------------------------------------------------------------------------------------------------------

# ==================== dynamic diagnostic plot =========================================================
# library(plotly)
plot_ly(ad100_qrs, x = ~ idx, y = ~signal, type = "scatter", linetype = 1, mode = "lines", showlegend = FALSE) %>% 
layout(xaxis = list(title = "sample (360 Hz)"),
yaxis = list(title = "unfiltered signal (mV)"),
legend = list(orientation = 'h')) %>%
add_markers(x = ~loc[!is.na(loc)], y = ~signal[!is.na(loc)], inherit = FALSE, name = "approximate location") %>% 
add_markers(x = ~loc_sr[!is.na(loc_sr)], y = ~signal[!is.na(loc_sr)], inherit = FALSE,
name = "corrected location")
# --------------------------------------------------------------------------------------------------------
