# Synchronized averaging
**Background:** Event related potentials (ERPs) are EEG signals that are obtained when the subject is asked to press a button when a certain shape appears on the screen (event).
This shape could appear in five regions: Upper right (UR), down right(DR), upper left (UL), down left (DL) and center of screen (CE). The characteristic template of the response of each region is buried in noise related to the overall background EEG activity of the brain.  
  
**Goal:** Removal of the EEG background noise to extract the charachterstic template of each response.  
  
**Approach:** Epochs of ERPs are obtained multiple times by repeated applications of stimulus; since the background noise is random with zero mean, they can be averaged out by using the onset of the stimulus to allign the epochs. This approach is known as synchronized averaging and it helps to increase the signal to noise ratio (SNR) when the noise and the desired signal's frequency spectrum overlap significantly. The performance of the algorithm is evaluated by measuring the SNR and the average Euclidean distance between the original nosie ERP signals and the averaged signals: 
  
**Resuls:**:

![alt text] (https://github.com/PhilippeMoussalli/Biomedical-data-analysis/blob/main/1)%20Synchronized%20averaging/figures/Extracted_temp.PNG "Title")


