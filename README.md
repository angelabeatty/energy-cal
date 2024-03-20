# energy spectrum calibration for .Chn and .spe data files
#energy_calibration.py is a python script that takes in .Chn and .spe spectrum data files and outputs the spectrum data in terms of energy. This script utilizes ortec_read read_chn script written by Oscar Tegmyr (https://github.com/tegmyr/ortec_read) to convert binary into a text file with spectrum data. \
The notebook contains examples analyzing spectrum data from sources of Am-241 and Cs-137, though other elements can also be used. \
\
#read_data takes in the data files and outputs the channel binning as x-data and spectrum signal as y-data. \
\
#gauss_opt takes in the y-data and upper/lower limit x-values corresponding to the range for the characteristic peak (note: currently this must be pre-determined by plotting the spectrum, future improvements would include a peak detection function), and applies scipy.optimize.curve_fit to fit a gaussian + step function to the peak. The output is the peak mean and the full-width half-max, both in terms of channel. \
\
#energy_opt takes in the channel peak, characteristic peak, and x-data for two sources. The channel and characteristic peaks are used to fit a linear function (energy = m*channel + b), the results being used to convert the x-data from channels to energies. The output is fit values and the x-data for both sources in terms of energy. \
\
#spect applies the gauss_opt function and has an optional plotting function so that the fit can be checked visually for accuracy. The output is the peak mean and the full-width half-max (FWHM), both in terms of channel.  \
\
#energy_fit takes in the outputs of #spect, the x- and y-data, and the characteristic peaks for two sources and applies the energy_opt function. The energy_opt function is applied, and the fitting results are used to recalculate the FWHM in terms of energy. The optional plotting function will plot both spectra in terms of energy and specify the characteristic peaks. The output is the fit values. 
