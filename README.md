# The-validation-of-our-proposed-method-of-the-gene-expression-binarization-using-Mathematica-

For the distribution to the community, the algorithm of our binarization method was implemented using both Mathematica and the R language (see, 
Gene-expression-binarization-R-programing).
Here, the Bi4Back package (the code of our binarization method using Mathematica) was implemented during my postdoctoral at IBISC laboratory in collaboration 
with M. Franck Delaplace. 

After, we tested our algorithm on real RNA-seq gene expression data using Mathematica and the R language. Here, we wanted to provide the codes (using Mathematica) for the validation of our proposed binarization method.  We first consider a simple example of an artificial
gene regulatory network that exhibits stable trajectories. This example is with a small dimension that permits first to present how we model a Boolean network with (ODE) systems. Then, artificial experiments data are generated using the simulation of the corresponding (ODE) system to be binarized using our proposed approach.  The outcomes are compared to binary profiles using the pre-set thresholds in the continuous models. Then, to test the performance of our proposed approach of binarization, we did the same thing considering well-known Boolean biological networks examples (with high dimension), taking Biane2018 Boolean network model, Sahin2009 Boolean network model, etc. Finally, we consider another simple example of 
an artificial gene regulatory network that exhibits stable orbits (oscillations) to show that our proposed approach of binarization is also working when the fluctuations exist. Using our binarization method, we have shown that the genes are correctly binarized based on the ODE simulations of artificial examples of gene regulatory networks or of well-known examples of Boolean biological networks.
