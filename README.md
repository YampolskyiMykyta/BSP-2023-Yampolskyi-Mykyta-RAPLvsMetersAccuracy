
# BSP-2023-Yampolskyi-Mykyta-RAPLvsMetersAccuracy  

The project aims to compare RAPL with Meter and access its accuracy.  
  
### Table of Contents  

- Introduction
- Installation  
- Usage  
- Contributing  
- License  
  
### Introduction  
Main project technical objectives:  
- Processing and organisation of data with the implementation of its visualisation.  
- Compare the energy consumption data obtained from RAPL and Meter measurements.  
- Perform statistical analysis to determine the level of accuracy and tolerance for RAPL measurements.  
- Have at the output a universal tool for comparing data with RAPL and Meter.  
  
### Installation  
To use this project, you will need to have Python, R and the required packages installed. Follow the steps below to set up the project:  
  
Clone the repository: git clone https://github.com/YampolskyiMykyta/BSP-2023-Yampolskyi-Mykyta-RAPLvsMetersAccuracy.git  
In Python:    
`pip install pandas`  
`pip install matplotlib` 
`pip install datetime`  
in R:    
Install the required packages by running install.packages(c("package1", "package2")), and run pre-build file.

### Usage

If you want to reproduce the resulting plots, the code is already in the repository, so you do not have to take the first step.  To run the code and reproduce data collection and analysis with your own Meter and RAPL values, follow these steps:

In Python for data collection:
 1. Load the required datasets, including the RAPL data, and Meter data.
 2. Call the `rewriteFileMETERandDeltas` function and providing the necessary arguments (once with FALSE and once with TRUE for two approaches). Then you will see how the data from the Meter will be parsed. 
 3. Call  the `get_everything_about_RMD` and `show_plot`functions to get information about the comparison between using conservative and less conservative approaches. And then create plots for Meter and RAPL data.

In Python for statistical analysis:
 1. Call the `get_speedup_rapl_median` function, providing the necessary arguments to see the speed-up test for the RAPL observations.
 
In R for statistical analysis:
 1. Run prebuild file to get all data in variables.
 2. If you want you can run descriptions and summaries for the data using `describe_orig_chunk` and `summary` functions.
 3. For normal distribution test use `get_some_nsw_data` function passing it the number of transactional models and threads.
 4. For correlation test use `get_pear_t_spear_f_plot_new` function passing it the number of transactional models.


### Contributing

Contributions to this project are welcome. If you encounter any issues or have suggestions for improvements, please submit a pull request or open an issue in the repository.

### License

This project is licensed under the [Apache2 license](http://www.apache.org/licenses/). Feel free to use and modify the code for your own purposes.