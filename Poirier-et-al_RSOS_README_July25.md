# Datafiles and Rcode for publication in Royal Society Open Science (July 2025)

[Access this dataset on Dryad](https://doi.org/10.5061/dryad.6t1g1jx9k)



## Description of the publication ##

Title: 
Urine washing and urinary odor profiles in relation to dominance rank status in wild male capuchin monkeys (<i>Cebus imitator</i>)

Authors: 
Alice C. Poirier(a,✭), Nelle K. Kulick(a,b), Suheidy Romero Morales(c), Marlen Kücklich(d,e), Brigitte M. Weiß(d,e), 
Claudia Birkemeyer(f,g,1), Anja Widdig(d,e,g,2), Amanda D. Melin(a,h,i,2), Katharine M. Jack(b,2,✭)
a Department of Anthropology and Archaeology, University of Calgary, Calgary, Canada
b Department of Anthropology, Tulane University, New Orleans, Louisiana, USA
c Área de Conservación Guanacaste, La Cruz, Costa Rica
d Behavioural Ecology Research Group, Faculty of Life Sciences, Institute of Biology, University of Leipzig, Leipzig, Germany 
e Department of Human Behavior, Ecology and Culture, Max Planck Institute for Evolutionary Anthropology, Leipzig, Germany
f Mass Spectrometry Research Group, Institute of Analytical Chemistry, University of Leipzig, Leipzig, Germany
g German Centre for Integrative Biodiversity Research (iDiv) Halle-Leipzig-Jena, 04103 Leipzig, Germany
h Department of Medical Genetics, University of Calgary, Calgary, Canada
i Alberta Children’s Hospital Research Institute, Calgary, Canada
1 Claudia Wiesner, prev. Birkemeyer
2 These authors contributed equally to this publication.
 

Abstract:
Urine plays an essential role in mammalian olfactory communication, though its potential role in primates has long been 
overlooked due to focus on their visual adaptations for communication. Here we combined behavioral and chemical data to 
test the role of urine in signaling male dominance in white-faced capuchins (Cebus imitator). We predicted that: 
1) urine washing (i.e., depositing urine onto hands/feet and rubbing them onto substrates) is more frequently performed by 
alpha than subordinate males, and 2) the chemical composition of alpha male urine is distinct from that of subordinates. 
We collected 457 hours of focal behavioral follows and 153 urine samples from 24 males in five groups at Sector Santa Rosa, 
Costa Rica. We extracted urinary volatile compounds into thermal desorption tubes, analyzed by gas chromatography-mass 
spectrometry. We found that alphas urine washed significantly more than subordinates, especially during the dry season when 
urinary odors can last longer and intergroup interactions are more frequent. We also found that dominance rank significantly 
predicted overall sample chemical dissimilarity. Our results support the hypothesis that urine may be an olfactory signaling 
medium; future experimental research is needed to test the extent to which urinary odors may be cues vs. evolved signals.

Keywords. Urine washing, olfactory communication, volatile organic compounds, chemosignaling, male dominance rank status, 
alternative male morphologies, primates.



## Description of the data and file structure ##

1. Poirier-et-al_RSOS_Dataset_FocalBehaviour_July2025 = Behavioral dataset
Variables: 
Date = focal data collection date; 
Focal_uniqueID = focal follow unique ID; 
Filename = focal data file per day and observer(sometimes including more than one follow);
Male = animal ID; 
DOB = date of birth; 
Rank = male rank; 
Group = monkey group ID; 
Observer = researcher who collected the data;
Focal_duration = duration of focal follow in H:MM:SS; 
UW_behav = UW event recorded (there can be more than one recorded during a follow; NA means no UW was recorded during the follow);
Seasonality = daily rainfall averaged over the 30 days prior to each date of our study period, continuous in mm;
Age = male age at the time of data collection, continuous in years.


2. Poirier-et-al_RSOS_Dataset_GCMS_SampleList_July2025 = sample list for chemical dataset
Variables:
Filename = urine sample unique ID;
Date = collection date;
Group = monkey group ID; 
Male = animal ID;
Rank = male rank;
Catch = catch type (direct or indirect);
Handler = researcher initials;
Batch = GC-MS batch;
Column = Column used on the instrument (old vs. new);
Age = male age at the time of data collection, continuous in years;
Volume = urine volume collected, continuous in mL;
SPG = urine mean specific gravity, continuous no unit;
Temp = ambient temperature at time of sample collection, continuous in C;
Seasonality = daily rainfall averaged over the 30 days prior to each date of our study period, continuous in mm.


3. Poirier-et-al_RSOS_Dataset_GCMS_CompoundList_July2025 = compound list for chemical dataset
Variables:
Comp = compound (peak) unique ID;
Filename = urine sample unique ID;
RT = peak retention time, continuous in min;
Area = peak area measured on the chemical profile, continuous no unit.

The sample list and compound list are combined by their common variable: Filename.


4. Poirier-et-al_RSOS_Rcode_July2025 = R code for statistical analyses
All the steps taken are described as comments in the code.



## Sharing/Access information ##

[Access this dataset on Dryad](https://doi.org/10.5061/dryad.6t1g1jx9k)

Link to other publicly accessible locations of the data:
[Github](https://github.com/AlicePoirier/Poirier-et-al_Urine-washing-and-urinary-odor-profiles-in-relation-to-dominance-rank-in-capuchins)

