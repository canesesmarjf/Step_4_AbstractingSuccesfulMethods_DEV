In this directory, we have the following subdirectories:
- bin/
- case_0/ (symlink)
- input_files/


The main workflow is as follows:
- Data produced by PICOS++ and found in case_0/ is used as the starting point
- That data is read in to MATLAB, plotted and them saved into a format that be passed to
  the C++ binary tree code (main_2.exe or main_3.exe). this data is saved in the input_files/
  directory.
- The c++ binary tree code is ran using ./run.sh 2 or 3 and the output it produces are called
  - ip_main_xx.csv
- This data is read back into MATLAB and contains indices to pick from the vectors x_p, v_p and a_p
- Data produced by C++ exectuable is plotted agains the PICOS++ data and compared

bin/ contains two binary files: (1) main_2.exe and (2) main_3.exe
What is the differenec between these two executables?
These files are actually symbolic links to the following locations
/home/jfcm/Repos/binaryTree/bin/main_2.exe
/home/jfcm/Repos/binaryTree/bin/main_3.exe

As it can be seen, these executables have been created in the following repo:

On branch multiDim_DEV
Your branch is up to date with 'origin/multiDim_DEV'.
origin	https://github.com/canesesmarjf/binaryTree.git (fetch)
origin	https://github.com/canesesmarjf/binaryTree.git (push)



case_0/ is a symbolic link to the PICOS++ results presented in ../Effort_1_Resampling/case_0
Where two species are used, one as a warm plasma while the other as a NBI species. This is the
directory used by the executables in bin/ in order to apply the  binary tree search.
