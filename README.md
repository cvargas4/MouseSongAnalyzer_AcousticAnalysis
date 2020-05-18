# MouseSongAnalyzer_AcousticAnalysis
Post-processing script for analyzing the Syllables.csv output from Mouse Song Analyzer used to analyzed mouse Ultrasonic Vocalizations (USVs).


## Set up and Running

The only dependancy is numpy

`
pip install numpy
`
<br/><br/>
<br/><br/>
To run the .py file, open terminal and enter:

`
python3 AcousticAnalysis_MSA.py
`
<br/><br/>
For the notebook, simply run the cell.
<br/><br/>
### Running

1. You will be prompted to enter a path for where the file you want to analyze resides.

`
  Folder/path of file to analyze: ~/MyFolder/MyCSV
`

2. Next, you will be prompted to input the name of the file.

`
  Name of .csv file to analyze: Syllables.csv
`
<br/><br/>
This will then take the .csv file, remove 'Unclassified' USVs, convert USVs with multiple jumps into 'M' type syllables (see table below for examples), and calculate means of acoustic features.
| Original Syllable Jumps| Converted Syllables| 
| -----------------------|:------------------:| 
|            s           |          s         |
|            d           |          d         |
|            u           |          u         |
|           ud           |          M         |
|           uu           |          M         |
|          udu           |          M         |

<br/><br/>
### Output
Script will create a folder named `MSAUSVAnalysis` with four .csv files.
*Features.csv*, *TransitionProbs_Overall.csv*, *TransitionProbs_Conditional.csv*
<br/><br/>
<br/><br/>
## Resources and References
Please see [MouseSongAnalyzer](https://github.com/cvargas4/MouseSongAnalyzer) to process .wav files of mouse USVs. This analysis will generate the file (ending in `Syllables.csv`) used by the script here.

This `AcousticAnalysis_MSA.py` python script replaces the`Acoustic_Analysis_Guided_v1.1.xlsx` file originally used in [Chabout et al., 2015](https://www.frontiersin.org/articles/10.3389/fnbeh.2015.00076/full) and [Chabout et al., 2016](https://www.frontiersin.org/articles/10.3389/fnbeh.2016.00197/full) for descriptive analyses of mouse USVs.

Please see those references for in-depth explanation of the transition and conditional probability calculations.
<br/><br/>
This python script was based on the above work by Jonathan Chabout, PhD, and written by Jannatul Ahmed (undergraduate student) and César Vargas (graduate student). 

This is research code. Please contact César for quesitons.
