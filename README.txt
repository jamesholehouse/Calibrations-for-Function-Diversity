The structure of this code is as follows. There are 4 folders in the highest level directory, which are,

1. Bacterial_Cells: containing the data analysis and calibration code and results for the bacterial cell protein abundance data.
2. Cities: containing the data analysis and calibration code and results for the city employee data.
3. Federal_Agencies: containing the data analysis and calibration code and results for the federal agency employee data.
4. Simulation_Code: containing the code required to simulate organizations from both naive and prescribed initial conditions.

The first 3 of these folders contains several elements:

a. A Julia file, e.g., calibration-cells.jl, that conducts the calibration procedure.
b. A Jupyter Notebook, e.g., Data-analysis-cells.ipynb, that shows how to access the data structures, and how to plot the diversity (number of unique elements) and the rank-frequency distributions.
c. A folder named "data" that contains the raw organization data.
d. A folder named "calibration_results" that contains the results of the calibrations.

Please note that these files will not work as downloaded, and the user will need to make sure that the correct locations in the filesystem are referenced at appropriate points in the code.

The calibrations were conducted on UNM's Center for Advanced Research Computing (CARC) supercomputers. Typically these calibrations were conducted in parallel.
To run this code, please make sure you have all the relevant Julia packages installed (I ran them using Julia 1.10.0) and use the following commands (respectively for 1, 2, 3 above):

1. julia -t 20 calibration-cells.jl 1 # where the "1" indicates the first data set in /data. It takes values 1, 2 or 3. 
2. julia -t 20 calibration-cities.jl
3. julia -t 20 calibration-federal-agencies.jl 1 #where the "1" indicates the first data set in /data. It takes values 1 or 2.

In the above code the "-t $NUM" refers to the number of cores that the code should be run on, since the main for loop conducting the calibrations is run in parallel.

Slurm scripts can be provided upon request. Please email at jamesholehouse1@gmail.com.
