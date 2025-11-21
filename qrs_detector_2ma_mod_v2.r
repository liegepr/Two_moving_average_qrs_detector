library(data.table)
library(gsignal)

# ============================= the main function ========================================================
ma_detector2 <- function(signal, srate = 360L, lowcut_f1 = 8L, highcut_f2 = 21L, filter_order = 3L, qrs_win1 = 35, 
    beat_win2 = 220, srate_ref = 360L, offset = 0.08, offset_win3 = 4L,
    slackness_red = FALSE, slackness_win1 = 0.200, slackness_win2 = 0.140,
    refractory_period = 0.3) {
    # signal: Numeric vector. ECG signal in mV.
    # srate: Integer. Sampling rate in Hz
    # lowcut_f1 & lowcut_f2: Integers. Bandpass for Butterworth filtering
    # filter_order: Integer. Default value for the Butterworth filter is set at 3th order (refer to Elgendi, 2013)
         # Default windows width as per Elgendi 2013 (optimized values)
    # qrs_win1: W1: Integer (nb. samples). "QRS window" ranges from 29 to 43 samples, for a sampling frequency (SF) of 360 Hz. Default value = 35.
    # beat_win2: W2: Integer (nb. samples). "Beat window", for a sampling frequency (SF) of 360 Hz, defaults to = 220.
        # In Porr & Powell 2019, window2 is multiple of window1 so that a window1 will not overlap two window2
        # In contrast, Elgendi, Jonkman & DeBoer 2010 & Elgendi 2013 allow for window1 overlapping two window2
    # srate_ref = Integer. Its value (360 Hz) should not be changed. 
         # 360 Hz is the frequency used in the three papers by Elgendi M. 
         # Remember that the optimized values for qrs_win1 and beat_win2 are valid for srate_ref = 360 Hz
    # offset: As scalar accounting for the threshold offset 'beta' (refer to Gradl & Elgendi 2015)
    # offset_win3: Integer (seconds). No optimized value. 
         # Gradl & Elgendi 2015 used 4s which was the maximum value for real-time calculations (embedded device)
    # slackness_red: Logical. slackness reduction as in Gradl & Elgendi 2015
    # slackness_win1: Scalar (seconds). Windows T1 for temporal correction (refer to Gradl & Elgendi 2015)
    # slackness_win2: Scalar (seconds). Windows T2 for temporal correction (refer to Gradl & Elgendi 2015)
    # slackness_win2: Scalar (seconds). Windows T2 for temporal correction (refer to Gradl & Elgendi 2015)
    # refractory_period: Scalar (seconds). Defaults to 300 ms (Porr & Powell 2019)

    nyquist_freq <- 0.5 * srate # Nyquist frequency
    low <- lowcut_f1 / nyquist_freq
    high <- highcut_f2 / nyquist_freq
  
    # signal filtering
        bandpass <- gsignal::butter(n = filter_order, w = c(low, high), type = "pass")
    signal_filt <- gsignal::filtfilt(bandpass, c(rep(signal[1],
        srate), signal, rep(signal[length(signal)], srate)))
    signal_filt <- signal_filt[(srate + 1):(length(signal_filt) -
        srate)]
    # signal squaring as in Elgendi 2013. Porr & Powell 2019 have been using abs()
    signal_squared <- signal_filt^2
    
    # Rolling means
    window1 <- trunc(qrs_win1 * srate / srate_ref)
        mwa_qrs <- frollmean(signal_squared, window1, align = "center") 
        # Center alignment as Elgendi 2013 (right alignment in Elgendi, Jonkman & DeBoer 2010 and Porr & Powell 2019)

    window2 <- trunc(beat_win2 * srate / srate_ref)
        mwa_beat <- frollmean(signal_squared, window2, align = "center")

     window3 <- trunc(offset_win3 * srate)
        mwa_noise <- frollmean(signal_squared, window3, align = "center")
    
    # If we stop there with running means, the 1st and last beats might not be identified
    # Because window2 is much larger than a qrs, too many running mean values are missing
    # The function below will substitute NA values with left-aligned and right-aligned  
    # running means at the left and right margins of the vector of means respectively.
    ma_fill <- function(mwa, signl, wind) {
        left_ma <- 1:(ceiling(wind / 2) - 1)
        right_ma <- (length(mwa) - floor(wind / 2) + 1):length(mwa)
        mwa[left_ma] <- frollmean(signl, wind, align = "left")[left_ma] 
        mwa[right_ma] <- frollmean(signl, wind, align = "right")[right_ma]
        mwa 
    }

    mwa_qrs <- ma_fill(mwa = mwa_qrs, signl = signal_squared, wind = window1)
    mwa_beat <- ma_fill(mwa = mwa_beat, signl = signal_squared, wind = window2)
    mwa_noise <- ma_fill(mwa = mwa_noise, signl = signal_squared, wind = window3)

    block <- fifelse(!is.na(mwa_beat) & mwa_qrs > (mwa_beat + mwa_noise * offset), 1, 0) 
    # the line above gives 0 if mwa_qrs <= mwa_beat
    # Now folding blocks and "silent" zones to get respective lengths
    
    blockt <- data.table(signal, signal_filt, signal_squared, mwa_qrs, mwa_beat, mwa_noise, block)
    
    # folding
    blocktc <- blockt[, .N, by=.(block, rleid(block))]
    
    # thresholding: any block smaller than window1 is excluded
    blocktc[block == 1L & N < window1, block := 0] 
    # block numbering
    blocktc[block == 1, block := frollsum(block, n=1:.N, adaptive = T) ]
    # zeroes are no longer needed; silent zones get NA
    blocktc[block == 0, block := NA]
    
    # unfold cleaned, numbered blocks
    blockt[, block.clean := blocktc[, rep(block, N)] ]
    # peak height on the squared scale
    blockt[, height := signal_squared]
    blockt[!is.na(block.clean), height := lapply(.SD,
                                                 FUN = function(x) fifelse(x == max(x), max(x), NA)),
           by = block.clean, .SDcols = "height"]
    blockt[is.na(block.clean), height := NA]
    
    # peak location
    blockt[, loc := 1:.N][is.na(height), loc := NA]
    
    # discard premature beats based on refractory period
    blockt[c(NA, diff(loc)) < refractory_period * srate, loc := NA]
    blockt[, idx := (1:.N)] # actual time index, used for merging data.tables, plotting, etc.

    # optional slackness reduction
    if (slackness_red) {
      srfun <- function(x, sig, t1, t2) {
        b1 <- as.integer(trunc(t1 / 2))  # we want t1/2 if t1 is an even integer, or its lower even integer if not
        b2 <- as.integer(trunc(t2 / 2))
        # the algorithm found Gradl & Elgendi 2015 is rather simple. However, the calculation will fail
        # at the 1st and last peaks if either the T1 or T2 window is not fully included 
        # in the sampling time interval
        # For T2, as straightforward solution was found
        # Our solution for T1 is more questionable, although it looks successful in most cases.
        T1 <- if ((x - b1) < 0) {
          c(1L, x + c(0L, b1))
        }
        else if ((x + b1) > length(sig)) {
          c(x + c(-b1, 0L), length(sig))
        }
        else {
          x + c(-b1, 0L, b1)
        }
        aR <- median(sig[T1], na.rm = TRUE)
        T2 <- x + (-b2):b2
        T2 <- T2[T2 >= 1L & T2 <= length(sig)] # discard T2 values outside from the range of samples
        # Gradl & Elgendi 2015 are using which.max(abs()), but this sometimes locate the S peak instead of the R one
        # Thus, we prefer using which.max()
        # yet to be checked: does omitting abs() work for all ECG leads?
        tU <- min(T2) + which.max(sig[T2] - aR) - 1L
        tU
      }
      swin1 <- slackness_win1 * srate
      swin2 <- slackness_win2 * srate

      loc2 <- blockt[!is.na(loc), .(sapply(loc, FUN = srfun, sig = blockt[, signal], t1 = swin1, t2 = swin2))]
      setnames(loc2, new = "idx")
      loc2[, loc_sr:=idx]
      # merging the corrected locations into the results table
      blockt <- loc2[blockt, on = c(idx = "idx")]
    }
    # time scale translation (if one needs to compare the results with those obtained using Matlab or C programs)
    # indexing starts from 1L in R, but from 0L in some other programs
    blockt[, idx:= idx - 1][, loc:= loc - 1]
    if ("loc_sr" %in% colnames(blockt)) {blockt[, loc_sr:= loc_sr - 1]}
    blockt[, time_stamp := idx / srate] 
    return(blockt)
}
# -----------------------------------------------------------------------------------------------------
ad100 <- fread("C:\\Users\\Philippe\\Documents\\Vélo\\Training\\Stats\\100.csv", sep = ",")
                      
# ==================== static diagnostic plots ==========================================================
ad100_qrs2 <- ma_detector2(signal = ad100[, MLII], srate = 360L, lowcut_f1 = 8L, highcut_f2 = 21L, filter_order = 3L, 
     qrs_win1 = 35L, beat_win2 = 220L, srate_ref = 360L, offset = 0.08, offset_win3 = 10L,
  slackness_red = TRUE, slackness_win1 = 0.2, slackness_win2 = 0.14, refractory_period = 0.3)

ad100_qrs20 <- ad100_qrs2[648600:650000]
ad100_qrs20[!is.na(block.clean), block.clean := 1]
ad100_qrs20[is.na(block.clean), block.clean := 0]

# plot the squared signal along with windows, blocks and estimated locations
plot(signal_squared ~ idx, data = ad100_qrs20, type = "l", xlab = "sample (360 Hz)")
lines(mwa_qrs ~ idx, data = ad100_qrs20, type = "l", col = 2) 
lines(mwa_beat ~ idx, data = ad100_qrs20, type = "l", col = 3) 
lines(I(0.2 * block.clean) ~ idx, data = ad100_qrs20, type = "l", col = 5) 
points(signal_squared ~ loc, data = ad100_qrs20, col = "darkred")
points(signal_squared ~ loc_sr, data = ad100_qrs20, col = "#036830", pch = 3)
legend(x = "topleft",          # Position
       title = "Slackness reduction",
       legend = c("No", "Yes", "W1", "W2", "Block"),  # Legend texts
       pch = c(1, 3, 0, 0, 0),           # Line types
       col = c("darkred", "#036830", 3, 2, 5),           # Line colors
       lty = c(0, 0, 1, 1, 1))                 # Line width

# same as above, using time as x-axis
plot(signal_squared ~ time_stamp, data = ad100_qrs20, type = "l", xlab = "time (s)")
lines(mwa_qrs ~ time_stamp, data = ad100_qrs20, type = "l", col = 2) 
lines(mwa_beat ~ time_stamp, data = ad100_qrs20, type = "l", col = 3) 
lines(I(0.2 * block.clean) ~ time_stamp, data = ad100_qrs20, type = "l", col = 5) 
points(signal_squared ~ I(loc / 360), data = ad100_qrs20, col = "darkred")
points(signal_squared ~ I(loc / 360), data = ad100_qrs20, col = "#036830", pch = 3)
legend(x = "topleft",          # Position
       title = "Slackness reduction",
       legend = c("No", "Yes", "W1", "W2", "Block"),  # Legend texts
       pch = c(1, 3, 0, 0, 0),           # Line types
       col = c("darkred", "#036830", 3, 2, 5),           # Line colors
       lty = c(0, 0, 1, 1, 1))                 # Line width

# plot the unfiltered signal and check location accuracy
# ad100_qrs20 <- ad100_qrs2[649730:649745]
# plot the unfiltered signal and check location accuracy
# ad100_qrs2[, loc_sr:=c(NA, loc_sr[1:649999])]
ad100_qrs21 <- ad100_qrs2[649500:650000]
plot(signal ~ idx, data = ad100_qrs21, type = "l", xlab = "sample (360 Hz)")
  points(signal ~ loc, data = ad100_qrs21, col = "darkred")
  points(signal ~ loc_sr, data = ad100_qrs21, col = "#036830", pch = 3)
  legend(x = "topleft",          # Position
       title = "Slackness reduction",
       legend = c("No", "Yes"),  # Legend texts
       pch = c(1, 3),           # Line types
       col = c("darkred", "#036830"),           # Line colors
       lty = 0)                 # Line width

# -------------------------------------------------------------------------------------------------------

# ==================== dynamic diagnostic plot =========================================================
library(plotly)
plot_ly(ad100_qrs20, x = ~ idx, y = ~signal, type = "scatter", linetype = 1, mode = "lines", showlegend = FALSE) %>% 
layout(xaxis = list(title = "sample (360 Hz)"),
yaxis = list(title = "unfiltered signal (mV)"),
legend = list(orientation = 'h')) %>%
add_markers(x = ~loc[!is.na(loc)], y = ~signal[!is.na(loc)], inherit = FALSE, name = "approximate location") %>% 
add_markers(x = ~loc_sr[!is.na(loc_sr)], y = ~signal[!is.na(loc_sr)], inherit = FALSE,
name = "corrected location")
# --------------------------------------------------------------------------------------------------------

