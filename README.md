# JWebbinar material

## Branch frozen for Webbinar 3

This Branch is a frozen version for JWebbinar 3 "JWST Pipeline in Imaging Mode". If you want a version of the materials that is being maintained and updates to the latest versions of the packages, see the `main` branch of this repository.

To download the material, you can clone this repository and navigate to this folder:  
`git clone https://github.com/spacetelescope/jwebbinar_prep.git ~/my_jwebbinar_prep`  
`cd ~/my_jwebbinar_prep`  
`git checkout webbinar3`  
`cd imaging_mode`  

Alternatively, you can manually download the individual notebooks by clicking on the imaging_mode folder, then click on each file, and right-click on "Raw" to select "Save link as".

To run these notebooks locally, you will need a Python environment set up that is capable of running Jupyter notebooks. If you are new to Python, you may want to look here: https://www.python.org/about/gettingstarted/ for some resources on how to install and learn the basics of Python. Depending on your preferences and system choices, you may find the install instructions there sufficient, but note that many scientists find it easier to use the Anaconda python distribution and package manager: https://www.anaconda.com/products/individual.

If you are using anaconda, you first create an environment with the name of your choice (we use `jwebbinar3` for this example):  
`conda create -n jwebbinar3 --file https://ssb.stsci.edu/releases/jwstdp/1.2.3/conda_python_macos-stable-deps.txt`  
`conda activate jwebbinar3`  
`pip install -r https://ssb.stsci.edu/releases/jwstdp/1.2.3/reqs_macos-stable-deps.txt`  
`pip install jupyter`  
(Note that there might be conflicts with the numpy version installed by the pipeline and by jupyter and you might have to downgrade jupyter.)

You can now run `jupyter notebook` in the folder where you have downloaded the webinar notebooks.

If you have problems or questions, feel free to reach out to the JWST Help Desk at http://jwsthelp.stsci.edu/.
