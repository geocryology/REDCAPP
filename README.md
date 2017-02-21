# REanalysis Downscaling Cold Air Pooling Parameterization (REDCAPP )

# What is REDCAPP
REDCAPP is a python-based open source software for parameterizing the temporal and spatial differentiation of surface effects and cold air pooling when downscaling reanalysis data in mountainous areas. REDCAPP can produce daily, high-resolution gridded fields of near-surface air temperature in mountains.

# Why REDCAPP is Powerful 
REDCAPP simulates high-resolution near-surface air temperaure well through providing a number of tools for

(1) downloading and manipulation of ERA-Interim data saved as netCDF4;

(2) contucting the interpoaltion of 2-meter air temperature and pressure level temperatures (or upper-air temperature);

(3) deriving proxy of land-surface effects from reanalysis data;

(4) addressing the spatially varying land-surface effects based on fine-scale DEM;


Additionally, the input data are not limited to ERA-Interim and could be extended to other reanalyses such as CFST, NCEP, MERRA or 20CRV2.

# How to Run REDCAPP
REDCAPP is wrotten by python (version 2.7) and public open. To run the software, please

(1) Get REDCAPP at https://github.com/geocryology/REDCAPP

(2) Make sure that the directory containing thie file (redcapp_example.py) is contained in your PYTHONPATH.

(3) Register to ECMWF (free) https://apps.ecmwf.int/registration/

(4) Follow the instructions for "Installing your API key" on
    https://software.ecmwf.int/wiki/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch

(5) Adapt the script below (settings and location) 

(6) Run the script

(7) Explore the results. Use a netcdf viewer to plot maps and time series.
    Panoply (https://www.giss.nasa.gov/tools/panoply) is a good one.

(8) Customise the code and use it for your project.



# Contact
Please let us know how things work. We hope this is useful for you.
Bin Cao (caobin198912@outlook.com)
