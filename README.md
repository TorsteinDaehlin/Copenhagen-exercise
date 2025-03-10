# Copenhagen exercise study
This repository contains code that computes segment and joint kinematics and kinetics from motion capture data. Specifically, this code is adapted to perform
inverse kinematics and dynamics on lower extremity motion capture data collected from participants performing the Copenhagen Adductor Exercise. However,
it would be relatively easy to extend this code for other tasks, marker sets, or outcomes of interest. 

If you use any part of the code shared in this repository for research, please cite: 

1. Dæhlin, T. E., Kennedy, M. D., Rouhani, H., & Chiu, L. Z. F. (2023). Effect of incline versus block heel-raise exercise on foot muscle strength and vertical jump performance – an 11-week randomized resistance training study. *Sports Biomechanics*, 1–17. [https://doi.org/10.1080/14763141.2023.2176353](https://doi.org/10.1080/14763141.2023.2176353)

### Required software
To run this code, you need [MATLAB](https://www.mathworks.com/products/matlab.html) and MATLAB's [Signal Processing Toolbox](https://www.mathworks.com/products/signal.html). The code herein, was originally developed in MATLAB 2019b, but has recently been heavily refactored for use in this research. The refactored code has only been tested in MATLAB 2024a and backcompatability cannot be guaranteed. Deprecated functions are functions from the original code base.

### How to set up and use the code
1. Fork this repository to your GitHub account.
2. Clone the fork to your machine. 

Once the code is cloned to your machine, you are ready to use the code to analyze your motion capture data.

### Structure of the code
The code processes the motion capture data in 5 main steps:
1. Essential preprocessing, including parsing and filtering marker and ground reaction force data, is performed in the first step. Code is included to parse input marker and force data collected using Qualisys Track Manager and exported using their `.mat` format. However, the code is structured such that further parser functions can be added seemlessly to the preprocessing pipeline.
2. The model is built based on a recorded standing trial which captures all anatomical and tracking markers for at least 1 frame. The `examples/model_setup.csv` file is used to specify which markers are attached to each segment. The model building subroutine defines a pelvis segment, thigh segment, leg segment, and foot segment. 
3. A 6 degree of freedom least square pose estimation algorithm is used to obtain segment positions and orientations. Segment angles and joint angles are calculated using Euler/Cardan angles. An XYZ rotation sequence was used in this study, corresponding to rotating about the medio-lateral axis, anterior-posterior axis, then vertical axis, as per the joint coordinate system convention.
4. Iterative Newton-Euler inverse dynamics is used to calculate net joint moments acting about each joint.
5. The data export function adds the outcomes to a `.mat` file which is reloaded at when starting the program. As such, the data can be processed all at one time, or over multiple sessions.

The code also contains functions that produces data-check figures, as well as some functions specific to the current project, such as an event definition function. In writing the code, I have attempted to plan with the easy of extending the code and allowing different uses in mind. 
