===============================================================================


# RobusTAD
A Tool for Robust Annotation of Topologically Associating Domain (TAD) Boundaries.

===============================================================================

### Objective:

RobusTAD calculates TAD boundary scores for every bin in the genome based on an interaction frequency matrix from Hi-C data. It also calls significant TAD boundaries.

##### Input: 

The only required input for RobusTAD is an Interaction frequency matrix for a given chromosome

__File properties:__


  * Tab-separated
  * Square matrix: number of rows should be equal to number of columns
  * Can contain raw or normalized counts
  * The file can have column and row names or not at all


An example input file is included in __Extras__.


__Example:__

_Option1_

```
	chr1-0	chr1-50 chr1-100	chr1-150	chr1-200
chr1-0	20	1	4	8	1	
chr1-50	5	25	8	2	0
chr1-100	4	7	18	3	6
chr1-150	3	2	7	30	8
chr1-200	5	2	1	9	27
```

_Option2_

```
20	1	4	8	1	
5	25	8	2	0
4	7	18	3	6
3	2	7	30	8
5	2	1	9	27
```




##### Output: 

RobusTAD outputs 2 files:

  * Boundary Scores: BoundaryScores_*.txt  
  * Significant Boundaries: TADBoundaryCalls_*.txt
  

BoundaryScores_*.txt  contains the Right, Left and Final scores for all bins in the provided IF matrix

__Example format:__

```
coordinates LeftBoundaryScore RightBoundaryScore TADscore
chr1-0 -0.27130659072985 -0.241287485278218 -0.241287485278218
chr1-50 -0.27130659072985 -0.241287485278218 -0.241287485278218
chr1-100 -0.27130659072985 -0.241287485278218 -0.241287485278218
chr1-150 -0.27130659072985 -0.241287485278218 -0.241287485278218

```

TADBoundaryCalls_*.txt  contains the Right, Left and Final scores for bins in the provided IF matrix that are called as TAD boundaries. Calls are made based on locating peaks in the boundary score profile that are above the set threshold.

__Example format:__

```
coordinates LeftBoundaryScore RightBoundaryScore TADscore
chr1-100 -0.27130659072985 -0.241287485278218 -0.241287485278218

```

Example output files are included in __Extras__.

##### Dependencies:

RobusTAD is written in R; a working R environement should be available.

RobusTAD also requires the "optparse" library to be able to parse command line options. You can install it in R using:


```
install.packages("optparse", repos="http://cran.us.r-project.org")

```



##### Usage: 
Rscript RobusTAD.R -i InputMatrix [options]


An example IF matrix is included in __Extras__. You can download it and test RobusTAD using:

```
Rscript RobusTAD.R -i Extras/IFmatrix_GM12878_Rao_Mbo_Chr20_50kb.txt
```

Example output files are also included in __Extras__.


##### Optional Parameters:

```

==============================================================================================================================================================

RobusTAD calculates TAD Boundary scores for each bin on a chromosome.
	Input: interaction frequency matrix for a chromosome.
	Output: 2 files: 
			I- file with TAD boundary scores (BoundaryScores_*): contains Right Boundary scores, Left Boundary scores and Final combined score (max(R, L)).
			II- file with TAD boundary calls (TADBoundaryCalls_*) identified by looking for local maxima above the set threshod.

Usage: RobusTAD.R -i InputMatrix [options]

==============================================================================================================================================================



Options:
	-i INPUT, --input=INPUT
		Interaction Frequency Matrix. Must be a square matrix: number of columns = number of rows

	-H, --header
		include -H if input contains a header/column names

	-n NORM, --norm=NORM
		indicates if IF matrix is raw or normalized [default = raw]; [options: {raw, norm}]

	-o OUTDIR, --outDir=OUTDIR
		output directory name

	-b BINSIZE, --binsize=BINSIZE
		binsize or resolution used in Hi-C analysis in kb [default = 50]

	-r MINRATIO, --minRatio=MINRATIO
		the ratio of TAD to background required to assign a good score [default = 1.5]

	-w MINWIN, --minWin=MINWIN
		minimum window around the bin used to calculate the TAD score in kb [default = 250]

	-W MAXWIN, --maxWin=MAXWIN
		maximum window around the bin used to calculate the TAD score in kb [default = 500]

	-T THRESHOLD, --threshold=THRESHOLD
		data percentile of TAD scores used to calculate threshold in order to call significant TAD boundaries. [default = 0.2]; [options: 0-1];
		the lower the threshold, the more stringent the TAD calls

	-h, --help
		Show this help message and exit



==============================================================================================================================================================

RobusTAD is available under a GPL liscence and comes with no warranties @ https://github.com/rdali/RobusTAD

==============================================================================================================================================================

```


##### License:

RobusTAD is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation.

RobusTAD is distributed in the hopes that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
