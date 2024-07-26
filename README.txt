Project Structure
This project analyzes the response data from different frequency experiments conducted on a structure. Below is the folder and file structure of the project:

1）Experiment_response_analysis - This folder contains the main script for analyzing the experiment response data.
2）Modal test frf - This folder contains the Frequency Response Function (FRF) data from modal testing.
Note: We provide only one optimal FRF because numerous experiments revealed that the dense frequency characteristics of the tensegrity structure often prevent obtaining a straightforward FRF. 
The specific parameter settings are configured through the BurstRandom module in the LMS system. 
Detailed parameters can be found in the Modal_test_frf text file. Running 'Stabilization.m' will generate the stabilization diagram.
3）Theoretical mode shape calculation - This folder contains the data and scripts related to the theoretical calculation of mode shapes.

Main Script
The main script for analyzing the experiment response data is located in the experiment_response_analysis folder and is named experiment_response_analysis.m.
Script Overview
The experiment_response_analysis.m script performs the following tasks:
Load Data: Loads the required data files, including the structural connectivity matrices, theoretical mode shapes, and experiment response data.
Select Experiment: Specifies which experiment's response data to analyze.
Process Data: Processes the selected experiment data to obtain the node displacements and reorders the node positions to match the connectivity matrices.
Filter Data: Applies a bandpass filter to the displacement data.
Analyze Frequency Response: Calculates the Fast Fourier Transform (FFT) of the filtered responses and plots the frequency response.
Calculate Mode Shapes: Computes the maximum amplitudes and phases of the node displacements.
Plot Results: Visualizes the experimental mode shapes and compares them to the theoretical mode shapes.

Mode_shape_calculation.m: This script is used to calculate the theoretical mode shapes of the structure. It includes:
Loading the structural connectivity matrices.
Defining material properties and geometric parameters.
Assembling the mass and stiffness matrices.
Performing eigenvalue analysis to obtain the natural frequencies and mode shapes.
Visualizing the mode shap

How to Use：
Open MATLAB.
Navigate to the experiment_response_analysis folder.
Open the experiment_response_analysis.m script.
Specify the desired experiment by changing the value of experiment_choice (1 for 7.28Hz, 2 for 9.7Hz, 3 for 10.45Hz) within the script.
Run the script.

Data Files
response_7_28hz.mat - Contains the response data for the 7.28Hz excitation.
response_9_7hz.mat - Contains the response data for the 9.7Hz excitation.
response_10_45hz.mat - Contains the response data for the 10.45Hz excitation.
Note: The original data was uploaded in an Excel format, which includes the parameters set through the motion capture system.

Functions
The function folder contains various utility functions used in the analysis. 
Ensure that this folder is added to the MATLAB path to enable the main script to access these functions.

Special Acknowledgment: 
Some functions are derived from the MOTES software published by Skelton and his team in 2019. 
We would like to express our special thanks to them! Readers can find related materials in the following paper:
Skelton, R., 2019. MOTES: Modeling of Tensegrity Structures. Journal of Open Source Software, 4(42), p.1613.


Contact
For any issues or questions, please contact the project maintainer.