# Copyright (C) 2023 Philippe Liège
# GPL GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

library(data.table)
library(gsignal)

# ============================= the main function ========================================================
ma_detector <- function(signal, srate = 360L, lowcut_f1 = 8L, highcut_f2 = 21L, filter_order = 3L, qrs_win1 = 35, 
    beat_win2 = 220, srate_ref = 360L, offset = 0.08, offset_win3 = 4L,
    slackness_red = FALSE, slackness_win1 = 0.200, slackness_win2 = 0.140,
    refractory_period = 0.3) {
    # signal: Numeric vector. ECG signal in mV.
    # srate: Integer. Sampling rate in Hz
    # lowcut_f1 & lowcut_f2: Integers. Bandpass for Butterworth filtering
    # filter_order: Integer. Default value for the Butterworth filter is set at 3th order (Elgendi, 2013)
         # Default windows width as per Elgendi 2013 (optimized values)
    # qrs_win1: W1: Integer (samples). 29 to 43 samples, for a sampling frequency (SF) of 360 Hz. Default = 35.
    # beat_win2: W2: Integer (samples). For a sampling frequency (SF) of 360 Hz. Default = 220.
        # In Porr & Powell 2019, window2 is multiple of window1 so that a window1 will not overlap two window2
        # In contrast, Elgendi, Jonkman & DeBoer 2010 & Elgendi 2013 allow for window1 overlapping two window2
    # srate_ref = Integer. Its value (360 Hz) should not be changed. 
         # 360 Hz is the frequency used in the three papers by Elgendi M. 
         # Remember that the optimized values for qrs_win1 and beat_win2 are valid for srate_ref = 360 Hz
    # offset:Scalar. Threshold offset beta (Gradl & Elgendi 2015)
    # offset_win3: Integer (seconds). No optimized value. 
         # Gradl & Elgendi 2015 used 4s which was the maximum value for real-time calculations (embedded device)
    # slackness_red: Logical. slackness reduction as in Gradl & Elgendi 2015
    # slackness_win1: Scalar (seconds). Windows T1 for temporal correction (Gradl & Elgendi 2015)
    # slackness_win2: Scalar (seconds). Windows T2 for temporal correction (Gradl & Elgendi 2015)
    # refractory_period: Scalar (seconds). Defaults to 300 ms (Porr & Powell 2019)

    nyquist_freq <- 0.5 * srate
    low <- lowcut_f1 / nyquist_freq
    high <- highcut_f2 / nyquist_freq
    # signal filtering
    #
    bandpass <- gsignal::butter(n = filter_order, w = c(low, high), type = "pass")
    signal_filt <- gsignal::filtfilt(bandpass, c(rep(signal[1],
        srate), signal, rep(signal[length(signal)], srate)))
    signal_filt <- signal_filt[(srate + 1):(length(signal_filt) -
        srate)]
    # signal squaring as in Elgendi 2013. Porr & Powell 2019 has been using abs()
    signal_squared <- signal_filt^2
    
    # Rolling means
    window1 <-  trunc(qrs_win1 * srate / srate_ref)
        mwa_qrs <- frollmean(signal_squared, window1, align = "center") 
        # Center alignment as Elgendi 2013 (right alignment in Elgendi, Jonkman & DeBoer 2010 and Porr & Powell 2019)

    window2 <-  trunc(beat_win2 * srate / srate_ref)
        mwa_beat <- frollmean(signal_squared, window2, align = "center")

     window3 <-  trunc(offset_win3 * srate)
        mwa_noise <- frollmean(signal_squared, window3, align = "center")
    
    # If we stop there with running means, the 1st and last beats won't be identified
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

    block <- data.table::fifelse(!is.na(mwa_beat) & mwa_qrs > (mwa_beat + mwa_noise * offset), 1, 0) 
    # the line above gives 0 if mwa_qrs <= mwa_beat
    # folding blocks and "silent" zones to get respective lengths
    block.le <- rle(block)
    # thresholding: any block smaller than window1 is excluded
    block.le$values[block.le$values == 1L & block.le$length < window1] <- 0
    # block numbering
    block.le$values[block.le$values == 1L] <- with(block.le, cumsum(values[values == 1L]))
    # zeroes are no longer needed; silent zones get NA
    block.le$values[block.le$values == 0L] <- NA
    # unfold cleaned, numbered blocks
    block.clean <- rep(block.le$values, block.le$lengths)
    # peak height on the squared scale
    height <- signal_squared
    height[!is.na(block.clean)] <- unlist(tapply(height, block.clean, 
        function(x) fifelse(x == max(x), max(x), as.numeric(NA))))
    height[is.na(block.clean)] <- NA
    # peak location
    loc <- which(!is.na(height))
    # discard premature beats based on refractory period
    loc[c(NA, diff(loc)) < refractory_period * srate] <- NA
    loc <- loc[!is.na(loc)]
    idx <- indx <- seq_along(height)
    indx[!indx %in% loc] <- NA
    loc <- indx

    if (slackness_red) {
    srfun <- function(x, sig, t1, t2) {
    b1 <- trunc(t1 / 2)  # we want t1/2 if t1 even, or its lower even integer if not
    b2 <- trunc(t2 / 2)
    T1 <- c(-b1, 0, b1)
    T2 <- (-b2):b2
    aR <- median(sig[x + T1], na.rm = TRUE)
    # Gradl & Elgendi 2015 are using which.max(abs()), but this sometimes locate the S peak instead of the R one
    # Thus, we prefer using which.max()
    # tU <- x - (b2 + 1) + which.max(abs(signal[x + T2] - aR))
    tU <- x - (b2 + 1) + which.max(signal[x + T2] - aR)
    tU
    }
    swin1 <- slackness_win1 * srate
    swin2 <- slackness_win2 * srate
        loc2 <- sapply(loc[!is.na(loc)], srfun, sig = signal, t1 = swin1, t2 = swin2)
    loc_sr <- rep(NA, length(loc))
    loc_sr[order(loc_sr) %in% loc2] <- loc2
    loc_sr <- as.integer(loc_sr)
    }

    dtb <- data.table(signal, signal_filt, signal_squared, mwa_qrs, mwa_beat, mwa_noise, block, block.clean, idx, loc, 
    height, loc_sr = if (exists("loc_sr")) loc_sr)
    dtb
}
# -----------------------------------------------------------------------------------------------------
ad100 <- fread("100.csv", sep = ",")
                      
# ==================== static diagnostic plot ==========================================================
ad100_qrs <- ma_detector(signal = ad100[, MLII], srate = 360L, lowcut_f1 = 8L, highcut_f2 = 21L, filter_order = 3L, 
     qrs_win1 = 35L, beat_win2 = 220L, srate_ref = 360L, offset = 0.08, offset_win3 = 10L,
    slackness_red = TRUE, slackness_win1 = 0.2, slackness_win2 = 0.14, refractory_period = 0.3)

ad100_qrs0 <- ad100_qrs[648000:650000]
ad100_qrs0[!is.na(block.clean), block.clean := 1]
ad100_qrs0[is.na(block.clean), block.clean := 0]

plot(signal_squared ~ idx, data = ad100_qrs0, type = "l", xlab = "sample (360 Hz)")
lines(mwa_qrs ~ idx, data = ad100_qrs0, type = "l", col = 2) 
lines(mwa_beat ~ idx, data = ad100_qrs0, type = "l", col = 3) 
lines(I(0.2 * block.clean) ~ idx, data = ad100_qrs0, type = "l", col = 5) 
points(signal_squared ~ loc, data = ad100_qrs0, col = "darkred")
points(signal_squared ~ loc_sr, data = ad100_qrs0, col = "#036830", pch = 3)
# -------------------------------------------------------------------------------------------------------

# ==================== dynamic diagnostic plot =========================================================
library(plotly)
plot_ly(ad100_qrs, x = ~idx, y = ~signal, type = "scatter", linetype = 1, mode = "lines", showlegend = FALSE) %>% 
layout(xaxis = list(title = "sample (360 Hz)"),
yaxis = list(title = "unfiltered signal (mV)"),
legend = list(orientation = 'h')) %>%
add_markers(x = ~loc[!is.na(loc)], y = ~signal[!is.na(loc)], inherit = FALSE, name = "approximate location") %>% 
add_markers(x = ~loc_sr[!is.na(loc_sr)], y = ~signal[!is.na(loc_sr)], inherit = FALSE,
name = "corrected location")
# --------------------------------------------------------------------------------------------------------

