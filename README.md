# registerTimeLapse
The batch script `registerTimeLapse.m` is designed to remove microscope jitter from a series of time lapse images.  Given a folder of multiple-image tif files, the script iteratively registers each image to the previous registered image.

## Dependencies
The function `dftregistration` is required for using this script and can be downloaded from the MATLAB file exchange:
> Manuel Guizar (2020). Efficient subpixel image registration by cross-correlation (https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation), MATLAB Central File Exchange. Retrieved December 8, 2020.

## User Inputs
The script requires the following user inputs:
- _imDir_ is the directory which contains the multiple-image tif files to be registered.
- _imTag_ is a string that identifies which images should be registered.  For example, this could specify a specific imaging channel.  To analyze all the images in the folder, set imTag to `''`.
- _saveDir_ is the location where the registered images will be saved.
- _overwrite_ allows the user to choose whether to overwrite or to skip previously registered image files.
- _usfac_ is the upsampling factor used to set the resolution of the registration (images are registered to within 1/usfac of a pixel).
- _otherChannels_ can be used to register other images with the same base name as the registered image.  If no other images need to be registered, set this option to `{}`.  Otherwise, specific a list of string that uniquely identify the other image sets to be image.  For example, if myImage_Ch1.tif is being registered using `imTag = 'Ch1'`, then myImage_Ch2.tif and myImage_Ch3.tif can be registered using the same offsets using `otherChannels = {'Ch2','Ch3'}`.
