# A two moving average qrs detector with optional temporal correction

* References

  * Elgendi, Mohamed; Eskofier, Björn; Dokos, Socrates; Abbott, Derek; Amaral, Luís A. Nunes. (2014): Revisiting QRS Detection Methodologies for Portable, Wearable, Battery-Operated, and Wireless ECG Systems. In PLoS ONE 9 (1), e84018. DOI: 10.1371/journal.pone.0084018.
  * Elgendi, Mohamed; Jonkman, Mirjam; Boer, Friso de (2009): Improved QRS Detection Algorithm using Dynamic Thresholds. In International Journal of Hybrid Information Technology 2, pp. 65–80.
  * Elgendi, Mohamed; Jonkman, Mirjam; DeBoer, Friso (Eds.) (2010): Frequency Bands Effects on QRS Detection. International Conference on Bio-inspired Systems and Signal Processing.
  * Elgendi, Mohamed; Talkachova, Alena (2013): Fast QRS Detection with an Optimized Knowledge-Based Method. Evaluation on 11 Standard ECG Databases. In PLoS ONE 8 (9), e73557. DOI: 10.1371/journal.pone.0073557.
  * Gradl, Stefan; Leutheuser, Heike; Elgendi, Mohamed; Lang, Nadine; Eskofier, Bjoern M. (2015): Temporal correction of detected R-peaks in ECG signals. A crucial step to improve QRS detection algorithms. In Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE Engineering in Medicine and Biology Society. Annual International Conference 2015, pp. 522–525. DOI: 10.1109/embc.2015.7318414.
  * Porr, Bernd; Howell, Luis (2019): R-peak detector stress test with a new noisy ECG database reveals significant performance differences amongst popular detectors: Cold Spring Harbor Laboratory.
  * Elgendi, Mohamed (2016): TERMA Framework for Biomedical Signal Analysis. An Economic-Inspired Approach. In Biosensors 6 (4). DOI: 10.3390/bios6040055.
  * Elgendi, Mohamed; Mohamed, Amr; Ward, Rabab (2017): Efficient ECG Compression and QRS Detection for E-Health Applications. In Sci Rep 7 (1), p. 459. DOI: 10.1038/s41598-017-00540-x.
  * Kristof, Florian; Kapsecker, Maximilian; Nissen, Leon; Brimicombe, James; Cowie, Martin R.; Ding, Zixuan et al. (2024): QRS detection in single-lead, telehealth electrocardiogram signals. Benchmarking open-source algorithms. In PLOS Digit Health 3 (8), e0000538. DOI: 10.1371/journal.pdig.0000538.
  * Porr, Bernd; Howell, Luis (2019): R-peak detector stress test with a new noisy ECG database reveals significant performance differences amongst popular detectors: Cold Spring Harbor Laboratory.
  * Porr, Bernd; Macfarlane, Peter W.; Naseer, Noman (2024): A new QRS detector stress test combining temporal jitter and F-score (JF) reveals significant performance differences amongst popular detectors. In PLoS ONE 19 (11), e0309739. DOI:      10.1371/journal.pone.0309739.
   * Reklewski, Wojciech; Miśkowicz, Marek; Augustyniak, Piotr (2024): QRS Detector Performance Evaluation Aware of Temporal Accuracy and Presence of Noise. In Sensors 24 (5), p. 1698. DOI: 10.3390/s24051698.

# Implementation details

* qrs_detector_2ma.r and qrs_detector_2ma_mod.r both return a data table containing the input signal, all intermediate calculations including filtered and squared signal, windows and blocks along with detected qrs (loc).
* One additional column (loc_sr) is reported if the user needs slackness reduction.
* The modified version (qrs_detector_2ma_mod.r) attempts to detect the very first and very last beats which are missed by design in the original version. The implementation of slackness reduction for these marginal peaks was found effective in most cases but though I am not saying it is 100% safe. 
* Two dependencies only: packages gsignal (butter() & filtfilt()) and data.table (fread(), data.table(), fifelse(), rleid(), frollingmeans() & foverlaps() + faster computation times) are required.
* No IF-THEN-ELSE, FOR or WHILE are being used.

# Validation against the Glasgow University Database (GUDB)

> High precision ECG Database with annotated R peaks, recorded and filmed under realistic conditions (https://researchdata.gla.ac.uk/716/ and https://github.com/berndporr/ECG-GUDB) 
> Howell, L. and Porr, B. (2018) High precision ECG Database with annotated R peaks, recorded and filmed under realistic conditions.
> University of Glascow
> Datacite DOI: 10.5525/gla.researchdata.716

* All recordings with annotations were tested: 123 for the chest strap and 106 for the loose cable setup.

* Without slackness correction
  * With zero tolerance

|   task    |  channel   |  TPR   |   F1   |
|:---------:|:----------:|:------:|:------:|
| hand_bike |   cable    | 0.3458 | 0.3454 |
| hand_bike | cheststrap | 0.3253 | 0.3236 |
|  jogging  |   cable    | 0.2661 | 0.2692 |
|  jogging  | cheststrap | 0.2797 | 0.2822 |
|   maths   |   cable    | 0.3487 | 0.3484 |
|   maths   | cheststrap | 0.4023 | 0.4022 |
|  sitting  |   cable    | 0.3394 | 0.3393 |
|  sitting  | cheststrap | 0.3922 | 0.3924 |
|  walking  |   cable    | 0.2908 | 0.2905 |
|  walking  | cheststrap | 0.3505 | 0.3504 |

  * Using a 40 ms tolerance window

> According to Porr & Howell, the default tolerance is a tenth of the sampling rate as may be read in the Physionet comparison algorithms. 
For a 250 Hz sampling rate, a tolerance of **40 ms** corresponds to (40e-3) * 250 = 10 samples. Porr & Howell also state that the widely used tolerance of **100 ms** which overestimates the performance of an algorithm should not be used. 
The detection is said to be right shifted (Porr & Howell 2019) but our sample plots show some cases of left-shifted detection. So, we used a tolerance interval which is symmetrical around the reference annotation.
The WFBD application guide (WAG.pdf) says that the match window specifies the maximum absolute difference in annotation times that is permitted for matching annotations. Its default value used by the bxb function is 0.15 seconds which is way too large.

|   task    |  channel   |  TPR   |   F1   |
|:---------:|:----------:|:------:|:------:|
| hand_bike |   cable    | 0.9949 | 0.9933 |
| hand_bike | cheststrap | 0.9846 | 0.9788 |
|  jogging  |   cable    | 0.9433 | 0.9519 |
|  jogging  | cheststrap | 0.9569 | 0.9628 |
|   maths   |   cable    | 0.9996 | 0.9988 |
|   maths   | cheststrap | 0.9665 | 0.9664 |
|  sitting  |   cable    |   1    | 0.9995 |
|  sitting  | cheststrap | 0.9678 | 0.9679 |
|  walking  |   cable    | 0.9869 | 0.9852 |
|  walking  | cheststrap | 0.9686 | 0.9681 |

* Using slackness correction

  * With zero tolerance

|   task    |  channel   |  TPR   |   F1   |
|:---------:|:----------:|:------:|:------:|
| hand_bike |   cable    | 0.9926 | 0.9911 |
| hand_bike | cheststrap | 0.9832 | 0.9774 |
|  jogging  |   cable    | 0.9381 | 0.9466 |
|  jogging  | cheststrap | 0.9529 | 0.9587 |
|   maths   |   cable    | 0.9994 | 0.9986 |
|   maths   | cheststrap | 0.9665 | 0.9664 |
|  sitting  |   cable    |   1    | 0.9995 |
|  sitting  | cheststrap | 0.9668 | 0.9669 |
|  walking  |   cable    | 0.9851 | 0.9835 |
|  walking  | cheststrap | 0.9678 | 0.9673 |

  * Using a 40 ms tolerance window

|   task    |  channel   |  TPR   |   F1   |
|:---------:|:----------:|:------:|:------:|
| hand_bike |   cable    | 0.9949 | 0.9933 |
| hand_bike | cheststrap | 0.9846 | 0.9788 |
|  jogging  |   cable    | 0.9433 | 0.9519 |
|  jogging  | cheststrap | 0.9569 | 0.9628 |
|   maths   |   cable    | 0.9996 | 0.9988 |
|   maths   | cheststrap | 0.9665 | 0.9664 |
|  sitting  |   cable    |   1    | 0.9995 |
|  sitting  | cheststrap | 0.9678 | 0.9679 |
|  walking  |   cable    | 0.9869 | 0.9852 |
|  walking  | cheststrap | 0.9686 | 0.9681 |


* Compare with the detect_rpeaks() function from the rsleep package

  * With zero tolerance

|   task    |  channel   |   TPR   |   F1    |
|:---------:|:----------:|:-------:|:-------:|
| hand_bike |   cable    | 0.1519  | 0.1516  |
| hand_bike | cheststrap | 0.2479  | 0.2444  |
|  jogging  |   cable    | 0.06779 | 0.07407 |
|  jogging  | cheststrap | 0.1936  | 0.1999  |
|   maths   |   cable    | 0.1692  | 0.1691  |
|   maths   | cheststrap | 0.3334  | 0.3272  |
|  sitting  |   cable    | 0.1922  | 0.1923  |
|  sitting  | cheststrap |  0.371  | 0.3538  |
|  walking  |   cable    | 0.07615 | 0.07591 |
|  walking  | cheststrap | 0.3008  | 0.2974  |

  * Using a 40 ms tolerance window
 
|   task    |  channel   |  TPR   |   F1   |
|:---------:|:----------:|:------:|:------:|
| hand_bike |   cable    | 0.9903 | 0.9905 |
| hand_bike | cheststrap | 0.9787 | 0.9552 |
|  jogging  |   cable    | 0.8355 | 0.8799 |
|  jogging  | cheststrap | 0.8662 | 0.9047 |
|   maths   |   cable    | 0.9989 | 0.9985 |
|   maths   | cheststrap | 0.9852 | 0.9662 |
|  sitting  |   cable    | 0.9981 | 0.9985 |
|  sitting  | cheststrap | 0.9952 | 0.9539 |
|  walking  |   cable    | 0.9897 |  0.99  |
|  walking  | cheststrap | 0.9933 | 0.9715 |

# Benchmarking
   * Time (ms) needed for processing 123 ECG recordings.

    > library(rbenchmark)
    > benchmark(
    >     dtbr <- apply(sample_list, 1, valid_fn, freq_sampling = 250L,  slred = TRUE)
    > )
    > test replications  elapsed  relative  user.self  sys.self  user.child  sys.child
    >             100   776.02         1     299.06     27.72          NA         NA
    

![plot1 *rightend*]

[plot1 *rightend*]: 2ma_detection_1.png "Windows, blocks & annotations"


![plot2 *detail*]

[plot2 *detail*]: 2ma_detection_2.png "Windows, blocks & annotations"

