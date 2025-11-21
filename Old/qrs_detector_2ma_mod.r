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
    indx <- seq_along(height)
    idx <- indx - 1L # actual time index, only used for plotting etc.
    indx[!indx %in% loc] <- NA
    loc <- indx

    if (slackness_red) {
    srfun <- function(x, sig, t1, t2) {
    b1 <- as.integer(trunc(t1 / 2))  # we want t1/2 if t1 even, or its lower even integer if not
    b2 <- as.integer(trunc(t2 / 2))
    # the algorithm found Gradl & Elgendi 2015 is very simple. However, the calculation will fail
    # at the 1st and last peaks if either the T1 or T2 window is not fully included 
    # in the sampling time interval
    # For T2, as straightforward solution was found
    # Our solution for T1 is more questionable, although is looks successful in most cases.
    T1 <- if ((x - b1) < 0) {
           c(1L, x + c(0L, b1))
            }
            else if ((x + b1) > length(sig)) {
            c(x + c(-b1, 0L), length(sig))
            }
            else {
            x + c(-b1, 0L, b1)
            }
    T2 <- x + (-b2):b2
    T2 <- T2[T2 >= 1L & T2 <= length(sig)] # discard T2 values outside from the range of samples
    aR <- median(sig[T1], na.rm = TRUE)
    # Gradl & Elgendi 2015 are using which.max(abs()), but this sometimes locate the S peak instead of the R one
    # Thus, we prefer using which.max()
    # yet to be checked: does omitting abs() work for all ECG leads?
    tU <- min(T2) + which.max(signal[T2] - aR) - 1L
    tU
    }
    swin1 <- slackness_win1 * srate
    swin2 <- slackness_win2 * srate
    loc2 <- sapply(loc[!is.na(loc)], srfun, sig = signal, t1 = swin1, t2 = swin2)
    loc_sr <- rep(NA, length(loc))
    loc_sr[order(loc_sr) %in% loc2] <- loc2
    loc <- loc - 1L # time scale translation (indexing starts from 1L in R)
    loc_sr <- loc_sr - 1L # time scale translation
    }
    time_stamp <- idx / srate
    dtb <- data.table(signal, signal_filt, signal_squared, mwa_qrs, mwa_beat, mwa_noise, block, block.clean, loc, idx,
    height, loc_sr = if (exists("loc_sr")) loc_sr, time_stamp)
    dtb
}
# -----------------------------------------------------------------------------------------------------

