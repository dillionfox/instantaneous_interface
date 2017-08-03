# instantaneous_interface
# python 2.7 (might work for other versions, haven't tried)
# Requirements: numpy, MDAnalysis

This code computes Willard-Chandler Instantaneous Interfaces as described in their paper, "Instantaneous Liquid Interfaces"

There are a few other codes out there that also accomplish this and are more efficient, but I found them difficult to use and not very flexible. This implementation is contained in one code and is designed to be readable and easy to modify. I have only tested this on a few systems, but it worked well for those. 

USEAGE: this code computes the interface around the heavy atoms of a given atom selection. You need to modify the variable called "selection_key" to identify your molecule(s) of interest according to MDAnalysis syntax. You also need to provide the path of the trajectory and structure files. I named the variables "DCD" and "PSF" but you can use any file types that work with MDAnalysis.

OUTPUT: The output is a pdb file for each frame in the trajectory containing the coordinates of the instantaneous interface.

If you don't need a very detailed surface and you want the code to run faster, try changing "dL" to 0.2, or even 0.5. This modification is not recommended for small molecules, however.

There are a few "joblib" functions included in the code that can parallelize some of the for loops which can be daunting when you have lots of heavy atoms, but in most cases they actually slow the code down so they are not used by default.

Feel free to share thoughts and comments.
