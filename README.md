#  Installation

Language: Python 3 (tested on 3.6, 3.8)

Direct Dependencies:

* numpy, pandas, matplotlib, seaborn, scipy, scikit-learn, ipython, jupyter
    
Special Dependencies:

* In-house library mesostat    
    * https://github.com/HelmchenLabSoftware/mesostat-dev
* Mesostat common dependencies
    * ipython, ipywidgets, h5py, json, numpy, scipy, scikit-learn, pandas, matplotlib, seaborn, pyqt5, networkx, statannot, patsy, shutil, opencv-python, multiprocessing, pathos, subprocess, psutil
* Mesostat special dependencies:
    * IDTxl:  https://github.com/pwollstadt/IDTxl
    * NPEET:  https://github.com/gregversteeg/NPEET

        
#  Structure

The front end of the code is contained in 6 jupyter notebooks inside the "src" folder

* "maria-analysis-1-explore.ipynb" - Visualization of raw data

* "maria-analysis-2-phases-and-intervals.ipynb" - Comparative analysis of task phases

* "maria-analysis-3-LR-CM.ipynb" - Statistical testing of differences between Left/Right turns and Correct/Mistake performance.

* "maria-analysis-4-trajectories-clustering.ipynb" - Study of neuronal trajectories in dimensionality-reduced manifolds, as well as activity clustering

* "maria-analysis-5-temporal-order.ipynb" - Study of within-period temporal orderability of neurons

* "maria-analysis-6-LR-CM-classification.ipynb" - Machine learning approach to classification of Left/Right turns and Correct/Mistake performance.



#  Usage
Each cell of each notebook computes a specific plot or a statistical result. Cells are labeled for clarity. The 2nd cell in each notebook contains paths to raw data (dff) and deconvolved data (deconv). These paths are absolute, and need to be adjusted by the user. Each path points to the root folder of the corresponding dataset.
