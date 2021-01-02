# Introduction
Coronaviruses (CoV) usually exhibit strong CpG deficiency on their genome. This coronavirus pattern of lacking CG dinucleotides can be explained by evasion of the immune response by the zinc finger antiviral protein (ZAP) in mammalian hosts. This genetic study tried to quantify the magnitude of CpG deficiencies, which can be served as an index of how potentially effective the coronavirus can escape from human’s immune system, among several coronavirus. They are severe acute respiratory syndrome coronavirus (SARS-CoV), severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2), middle east respiratory syndrome coronavirus (MERS-CoV) and two bat-derived severe acute respiratory syndrome-like coronaviruses (bat-SL-CoVZC45 and bat-SL-CoVZXC21).

# Methods
By using the Biostrings package, I calculated the occurrences and percentages of C, G and CG dinucleotides of each coronavirus whole-genome sequence. Then I split each sequence into 29 to 30 1000 bp windows based on their sequence length. I calculated both the GC content (pC + pG, where pC and pG are the percentages of bases being G or C) and the CpG deficiency (observed-to-expected CG ratio: calculated by pCG/(pC*pG), pCG is the percent of CG dinucleotides occurrence) on each 1000 bp window for each coronavirus sequence, and generated the density plots and scatter plot of GC contents and CpG deficiency for all these 5 virus sequences to see if there is any heterogeneity of CpG deficiency among these coronaviruses. 

# Main Result
![alt text](https://github.com/Holin-Chen/coronavirus-CpG-genetic-analysis/blob/main/plots/I_CpG%20plot.png)

The curves of GC contents and CpG deficiency for SAR-CoV-2 reach peaks the most quickly than other 4 coronaviruses, meaning the CpG deficiency for SAR-CoV-2 is the most apparent among these 5 coronaviruses. The curve of CpG deficiency for MERS-CoV reaches peak at last in high density, indicating the CpG deficiency for MERS-CoV is relatively moderate comparing to other coronaviruses.

![alt text](https://github.com/Holin-Chen/coronavirus-CpG-genetic-analysis/blob/main/plots/scatter%20plot.png)

most scatters of SAR-CoV-2 with its mean point (pCG = 0.3786, CpG deficiency=0.3803) locate at the bottom left area of the plot, proving again that SAR-CoV-2 has the lowest CG content and strongest CpG deficiency among the five coronaviruses. The mean CpG deficiencies among 1000 bp windows for SAR-CoV, MERS-CoV, bat-SL-CoVZC45 and bat-SL-CoVZXC21 are 0.4432, 0.5577, 0.4224, 0.4224 separately, which are all larger than that in SAR-CoV-2 of 0.3803. 
