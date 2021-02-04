# Removal of powerline noise interfernce 

## Background:
A major source of interference in ECG signals is the 50 or 60 Hz power-line frequency. The power-line frequency component can be removed by using a notch filter, in which a zero is placed on the unit circle at the location corresponding to this frequency.
h<sub>&theta;</sub>(x) = &theta;<sub>o</sub> x + &theta;<sub>1</sub>x
  
## Goal:
Removal of the EEG background noise to extract the characteristic template of each response.  
  
## Approach:
Epochs of ERPs are obtained multiple times by repeated applications of stimulus; since the background noise is random with zero mean, they can be averaged out by using the onset of the stimulus to allign the epochs. This approach is known as synchronized averaging and it helps to increase the signal to noise ratio (SNR) when the noise and the desired signal's frequency spectrum overlap significantly. The performance of the algorithm is evaluated by measuring the SNR and the average Euclidean distance between the original nosie ERP signals and the averaged signals: 
  
## Resuls:
  
**Extracted template overalapped on all ERP epochs**
![image info](./figures/temp_overlap_ERP.PNG)  
  
**Extracted templates visualized individually**
![image info](./figures/Extracted_temp.PNG)

**SNR and Euclidean distance**  
  
![image info](./figures/results.PNG)




