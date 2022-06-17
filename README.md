# JWebbinar material

## Branch frozen for Webbinar 4

This Branch is a frozen version for JWebbinar 4 "Pipeline in spectroscopic mode". If you want a version of the materials that is being maintained and updated to the latest version of the packages, see the `main` branch of this repository.

To download the material, you can clone this repository and navigate to this folder:  
`git clone https://github.com/spacetelescope/jwebbinar_prep.git ~/my_jwebbinar_prep`  
`cd ~/my_jwebbinar_prep`  
`git checkout webbinar4`  
`cd spec_mode`  

Alternatively, you can manually download the individual notebooks by clicking on the spec_mode folder, then click on each file, and right-click on "Raw" to select "Save link as".

To run these notebooks locally, you will need a Python environment set up that is capable of running Jupyter notebooks. If you are new to Python, you may want to look here: https://www.python.org/about/gettingstarted/ for some resources on how to install and learn the basics of Python. Depending on your preferences and system choices, you may find the install instructions there sufficient, but note that many scientists find it easier to use the Anaconda python distribution and package manager: https://www.anaconda.com/products/individual.

Using anaconda o miniconda, you first create an environment with the name of your choice (we use jwebbinar4 for this example):
MacOS:
`conda create -n jwebbinar4 --file https://ssb.stsci.edu/releases/jwstdp/1.1.0/conda_python_macos-stable-deps.txt`  
`conda activate jwebbinar4`  
`pip install -r https://ssb.stsci.edu/releases/jwstdp/1.1.0/reqs_macos-stable-deps.txt`  
`pip install jupyter`  
Linux:  
`conda create -n jwebbinar4 --file https://ssb.stsci.edu/releases/jwstdp/1.1.0/conda_python_stable-deps.txt`  
`conda activate jwebbinar4`  
`pip install -r https://ssb.stsci.edu/releases/jwstdp/1.1.0/reqs_stable-deps.txt`  
`pip install jupyter`  
(Note that the there might be a conflict with the version of numpy and you might have to downgrade jupyter.)

You can now run `jupyter notebook` in the folder where you have downloaded the webinar notebooks.

If you have problems or questions, feel free to reach out to the JWST Help Desk at http://jwsthelp.stsci.edu/.
