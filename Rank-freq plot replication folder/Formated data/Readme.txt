Fed_rank_freq.csv: 
Rank-freq of all agencies bundled together. 

"OCC" occupation category code
"n" number of people in this occupation category 
"rank" the rank of occupation 
(ignore first col, which is python DF index)


fed_rank_freq_each_org.csv:
Rank-freq of each agency:
 "AGY": code for federal agency (DOD, DOJ,... )
"OCC": code for occupation category. 
"EMPLOYMENT": number of people in this category. 


msa_rank_freq_each_city.csv:
rank_freq_for each city
(ignore first col, which is python DF index)

AREA: code for MSA (Metro Statistical Area)
OCC_CODE: occupation code
AREA_NAME: name of Metro area
TOT_EMP: number of people employed in the s occupation code. 

uni_rank_freq_each_uni.csv:
Universities offering Bachelor's degrees and above. 
DegreeID: Identifiere of the degree. First number represents the major. The second represents the level of degree (Bachelor's, Master's, PhDs. )
UNITID: ID for universities
nDegrees: number of degreed awarded in that degree in a year. 


uni_assoc_rank_freq_each_uni.csv
Same thing as above but for associate level universities
There's one degree level in this. 

Norway.csv
N: Total employees of company
D: Number of distinct occupations of company. 
This is the aggregated data of 5 companies. 

Fed.csv:
Diversity scaling for Fed agencies. Col names same as Norway. 

Uni.csv
Diversity scaling for Bachelor + universities

uni_associate.csv
Diversity scaling for associate schools. 

msa.csv
Diversity scaling for cities. 



