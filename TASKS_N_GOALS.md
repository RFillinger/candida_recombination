# To Do's:

- # 2
- # 4

- # 10

- # 19
- # 20
- # 21
- # 22
- # 23
- # 24

- # 30 
- # 31
- # 32
- # 33

1. ~~Recombination along the chromosome relative to the centromere.~~  ~~**Still need to do the chr5 recombination part of this figure.**~~ **Vetted all the markers in centromeres; no recombinatons in them.**

2. Make the mitochondrial main figure panel circular. 

3. ~~How many SNPs between SC5314 and P60002 and SC5314 and 529L~~  **157,900 (yes, that is the exact number) unique SNPs  and INDELS between MAY1 and 529L, not including the mitochondria. 73,518 unique SNPs and INDELs between MAY1 and P60002, also not including the mitochondria**

4. Need images of opaques for mating parents
    1. ~~Plate strains from storage to YPD (MAY 73, MAY 85)~~
    2. ~~Streak onto SCD~~
    3. Photograph opaque sectors
    4. Photograph opaque cells under microscope

5. ~~Calculate similarity in aneuploid profile for each mating pair (Figure 1C)~~ **No difference between distributions of aneuploidy (Kruskal-Wallis:  chi-squared = 0.0116099, p-val = 0.91419)**

6. ~~Identify the range of SC5314, 529L, and heterozygote marker position in SC534x529L progeny~~ **( 529L x SC: HET MAX: 99.5%, MIN: 0%, P1 MAX: 100%, MIN: 0%, P2 MAX: 40.4%, MIN: 0%  ) (P60002 x SC: HET MAX: 2.6%, MIN: 0.216%, P1 MAX: 99.8%, MIN: 97.4%, P2 MAX: 0.43%, MIN: 0% )**

7. ~~Determine the average percent het and P60002 in the SCxP6 progeny~~ **(529LxSC ddRAD-Seq: "1"-4524, "2"-1521, "n"-23666, total_informative-29711, 79.6% Het; 
P6xSC ddRAD-Seq: "1"-30821, "2"- 20, "n"-395, total_informative-31236, 1.3% Het)
529LxSC WGS: "1"-10759, "2"-8544, "n"-73563, total_informative-92866, 79.2% Het**

8. ~~Number of recombination events in each of the two matings (upper and lower bounds)~~ **Max: 102 recombinations, Min: 0 recombinations**

9. ~~Is the distribution of recombination tract lengths a Poisson for either mating?~~ **This one is more complicated than initially thought. It turns out pretty much any unimodal distribution is a Poisson distribution and it's a matter of calculating the lambda coefficient. The smaller the lambda, the more likely that smaller values will be more frequent on the distribution. SCx529L lambda = 0.506, SCxP60002 lambda = 0.506; both of these distributions of recombinant fragment sizes are front-heavy and likely to be smaller, which makes sense. This was calculated by graphing log(fragment size frequency) by log(fragment size) and fitting a linear model. The slope of the linear model is lambda. I had to compress the data down to just 10 bins because log functions don't play well with zeroes, so I'm unsure of the quality of the lambdas that I measured, but it jives with what we've got.**
    1. ddRAD: SCx529L lambda = 0.506
    2. ddRAD: SCxP60002 lambda = 0.506
    3. WGS:   SCx529L lambda = 0.509

10. PCR and sanger sequencing of Mitochondria for additional evidence of recombination
    1. ~~Design primers of 3 regions in recombinant mitochondria~~
    2. PCR regions out
    3. Submit for sequencing
    4. Integrate sequences into a supplemental figure

11. ~~Ensure compared frequencies of aneuploids is based on fraction/percent of total progeny, not absolute numbers. Need W stat for KW test of equal aneuploid frequency.~~ **Aneuploid frequencies were normalized to progeny count W: 0.0116, P.value: 0.9141**

12. ~~Number of recombination events in the SCxP6 mating (upper and lower bounds)~~ **(Max: 31, Min: 0, Mean: 3.91)**

13. ~~What is the maximum tract length for each of the two matings~~ **(SCx529L: max: 1604219bp, min: 3bp SCxP60002: max: 1582704bp, min: 3bp ), 3bp is the minimum distance between two markers, so the actual is probably larger, but it's occurred between two markers that are registered 3bp apart**

14. ~~Two-way ANOVA comparing the two parasexual progeny groups for % heterozygous and homozygous determined for each progeny genome)~~ **I did two one-way anovas looking at heterozygous and homozygous counts separated by matings. Both were extremely significant between SCx529L and SCxP6 progeny: HETEROZYGOUS f-statistic = 316.4, p.value < 2e-16, HOMOZYGOUS f-statistic = 378.2, p.value < 2e-16**

15. ~~Summed relative frequency of recombination divided by the number of bins for each chromosome for WGS analysis~~ **Done! It's in a supplemental file titled Recombination_Density_Per_Chromosome**

16. ~~What is the range of and the average number of recombination events in the 5 WGS progeny?~~ **Max: 505, Min: 430, Mean: 458.8**

17. ~~What is the  average number of recombination events in the SCx529L ddRAD-Seq progeny?~~ **Max: 102, Min: 1, Mean: 35.04**

18. ~~What is the range of recombination tract lengths in the WGS data?~~ **Max size: 724659bp, Min size: 18bp**

19. Update Figure 1.

20. Update Figure 2.

21. Update Figure 3.

22. Update Figure 4.

23. Compile supplemental figures with their main-figure counterparts

24. Write methods section

25. ~~Inform Anna and Abhishek of their respective method section writing duties~~

26. ~~We need to know the number of progeny for each mating (SCx529L and SCxP6).~~ **SCx529L: 75 Progeny, SCxP60002: 70 Progeny**

27. ~~What is the tract length for the WGS dataset for SCx529L progeny (min/max).~~ **MAX: 724659 bp, MIN: 1 bp**

28. ~~Check a few of the 1+ nt recombinations in the WGS data.~~ **This is a long story... but a complete one! I rewrote the pipeline to deal with strange markers and replacement markers that have been manually curated** 

29. ~~The number of markers for ddRAD-Seq of both matings and for WGS analysis in SCx529L~~ **(SCxP6 ddRAD: 740, SCx529L ddRAD: 678, SCx529L WGS: 18623)**

30. Break the variants between SC5314 and 529L or P6002 into SNPs and indels. How many of each?

31. LOH recombination figure 

32. Rewrite elements of the paper to include this data. 

33. Graph in Figure 4 is missing data; fix this in the `phenotypes_graph_maker.R` program

