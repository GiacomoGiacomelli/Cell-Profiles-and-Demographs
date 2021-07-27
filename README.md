# Cell profiles and demographs

## Available Files:
ProfilingCells_EPI.ijm
- Input: Multi Channel Microscopy Image (Any Image)
- Input: RoiSet.zip (Segmented Lines with defined Width -> User Defined)
- Output1: "Gray*.txt" (Position X, Channel 1 Intensity from ROI #*) 
- Output2: "Red*.txt" (Position X, Channel 2 Intensity from ROI #*) 
- Output3: "Blue*.txt" (Position X, Channel 3 Intensity from ROI #*)
- Extra requirements 1: Create a new folder named "Profiles" in the same folder containing the image
- Extra requirements 2: Rename the folder after completing an image "Profiles" -> "ProfilesN" (Repeat "Extra requirements 1 if analysing a second image")

CPAD1_Extract_Cell_length_from_membrane_stain_V1.R
- Input: "Profiles*" folders containg fluorescence profiles of membrane stained cells (Usually "Red*.txt")
- Ouput: Curated Profiles ("Gray*.txt","Red*.txt","Blue*.txt")
- Extra requirements: Create a "Membranes*" with the same numbering of the "Profiles*" folders

CPAD2_Cell_Profiles_And_Demographs_V3.R
NOTE: This Script uses the "findpeaks" function from -> https://rdrr.io/bioc/alsace/man/fitpeaks.html
- Input1: "Membranes*" or "Profiles*" folders containing fluorescence profiles ("Gray*.txt", "Red*.txt, "Blue*.txt").
- User Defined1: PlotType ("CellNorm/PopNorm") -> determines the type of normalization
- User Defined2: m*100/mayim OR 1+(m_c*100/mayim_c)  -> determines the type of normalization
- Output1: "../Profiles_ordered.txt" -> file containing the fluorescence profiles ordered by cell length
- Output2: "../Profiles_ordered_matrix.txt" -> matrix used to represent the fluorescence profiles as demographs via hist2d()
- Output3: "../maxim.txt" -> maximum length among the profiles
- Output4: "../Demograph.png" -> demograph
- Output5: "../Cell*.png" -> Fluorescence profiles and peak position
- Output6: "../CellLength_and_Peaks.txt" ->  file containing cell length and the number of peaks contained in each fluorescence channel 
