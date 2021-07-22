# Monkilock, Wu et al., 2021
This repository contains the codes and data for the Monkilock Project.

### Codes:
- **Monkilock_Data_Processing.ipynb**: In this notebook, we processed raw fixation data of 5 macaques and explored the relationship between their fixation patterns and the "surprisal" of events in a trial
- **Analysis_Code.Rmd**: Data processing, modeling, visualization
- **helper-lib.R**: helper functions used in Analysis_Code.Rmd

### Data:
- **cdv-timestamp folder**: one file for each trial; each file contains one timestamp per line
- **csv folder**: one file for each trial; each file contains one event (pop-up) per line
- **csv-combined.csv**: one event per line for all trial files (aggregated csv file)
- **csv-surprisal-prob.csv**: output file with surprisal value generated by Monkilock_Data_Processing.ipynb
- **csv-surprisal-prob_updated.csv**: ready for analysis data file generated by Analysis_Code.Rmd
- **AllSeq.csv**: include event information of all 80 sequences






