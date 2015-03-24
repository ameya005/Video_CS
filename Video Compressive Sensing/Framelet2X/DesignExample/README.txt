

These MATLAB programs reproduce the design in Example 1
of the paper:

 I. W. Selesnick and A. Farras Abdelnour,
 Symmetric wavelet tight frames with two generators,
 Applied and Computational Harmonic Analalysis, 
 17(2):211-225, September 2004.
 (Special Issue on Frames in Harmonic Analysis, Part II).

Function list (see individual functions for more detail):

DesignExample - This is the main program that executes all
the steps in the design procedure. The remaining programs
are 'helper' functions.

binom - Binomial coefficients. 

dbleroots - Finds the roots of a polynomial having all roots
of even degree. This function is more numerically accurate
than the 'roots' command in MATLAB because it takes the
derivative first.

extractf - Given H(z) and P(z), this function finds a polynomial
F(z) such that H(z) = F(z) P(z) if such a function exists.
This function is more numerically accurate than the 'deconv'
function in MATLAB.

flip - Reverses the order of elements in a vector.

scalfn - Computes samples of a scaling function.

up - Upsamples a vector (inserts zeros between each element)

wletfn - Computes samples of a wavelet.

z2x - Implements the change of variables x = (-z + 2 - 1/z)/4
that converts an odd-length symmetric filter into a generic
polynomial.





