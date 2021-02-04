# Independant Component Analysis (ICA) 

## Background:
Independent component analysis (ICA) is a signal processing tool used to separate an instantaneous mix of independent, non-Gaussian random variables. This is used mainly when we have different sensors that are picking up different signals and we want to decompose the signal into its original components to uncover the different sources or to denoise the original mixtures. The underlying assumption is that the sources are linearly mixed, do not have a Gaussian distribution and are independent one from another. The mixture can be represented with the following mixing model:
**X=AS**  
  where:
  *	**X:** is the mixture, that is the signal that is captured. Its rows are equal to the number of sensors that the signal is captured from.  
  * **A:** is the mixing matrix. Each column of the mixing matrix contains the mixing weights for a defined signal source.
  * **S:** S is the source matrix. Each row of this source matrix contains an independent component that is to be estimated. The number of components that can be derived cannot exceed the number of available sensors.
    
 Thus, the estimation of the sources can be simplified by the following equation: 
 **Y=WX** 
 
 Where Y is our signal estimation (Y≈S). and W is the inverse or pseudo-inverse of the mixing matrix A also known as the de-mixing matrix (W≈A<sup>-1</sup>). In other words, if we know the linear combination by which the sources have been mixed, we can recover the original sources present in the signal. The steps by which it achieves that are summarized as follows:  
   
 **1) Whitening:** Before separating the independent components, the mixture X is whitened so that its different components are made uncorrelated to each other. Its aim is to recover the original shape of the data that was altered during the mixing process. This is usually done with projecting the data unto its principal components as to obtain components that are uncorrelated with one with the other (covariance matrix of unity). Whitening the data renders the estimation of the independent components much easier to derive.  
   
   **2) FastICA algorithm:** After obtaining the uncorrelated sources. We are left with the task of reducing the Gaussianity of the data by maximizing the negentropy of the whitened data. In other words, we want to obtain the estimated sources Y that are statistically independent from each other. This follows the central limit theorem dictates that any mixture of random variables is more Gaussian that the original variables. At each iteration of the algorithm, the values of W are updated as to make the sources more independent from each other (less Gaussian). When W converges, we obtain the estimated independent source.

  
## Goal:
explore different ICA properties by mixing four different mixtures with different types of mixtures and noises and try to denoise specific sources from our mixture. The aim is to get acquainted with the way ICA works and specifically the FastICA algorithm in that context.
  
## Approach:
The power-line frequency component can be removed by using a notch filter, in which a zero is placed on the unit circle at the location corresponding to this frequency. The notch filter can be derived from the angular position of the required zero on the unit cirlce as given by:  
  
  
## Results:
  
**Power Spectral density of ECG signal before and after filtering**

  
**Filtered signal overlapped over the original signal**







