Procedure:

1. download the experimental branch of iCIPT2
2. compile iCIPT2
3. set the env 
   e.g. 
export ICI_CPP        =yourpath/iCIPT2_CXX/bin/iCI_CPP_NEW.exe
export ICI_CSF_CPP    =yourpath/iCIPT2_CXX/bin/iCI_CPP_NEW.exe
export ICI_DET_CPP    =yourpath/iCIPT2_CXX/bin/iCIPT2_D2h_Det.exe
export ICI_DET_COOV   =yourpath/iCIPT2_CXX/bin/iCIPT2_Coov_Det.exe
export ICI_DET_DOOH   =yourpath/iCIPT2_CXX/bin/iCIPT2_Dooh_Det.exe
export ICI_CSF_COOV   =yourpath/iCIPT2_CXX/bin/iCIPT2_Coov_CSF.exe
export ICI_CSF_DOOH   =yourpath/iCIPT2_CXX/bin/iCIPT2_Dooh_CSF.exe
export ICI_CSF_CVS_CPP=yourpath/iCIPT2_CXX/bin/iCIPT2_D2h_CSF_CoreExt.exe
4. run 00-integrals.py to generate  the integrals.

python3 00-integrals.py --basis=BASIS

In this test, BASIS should be one of ccpvdz ccpvtz ccpvqz ccpv5z

NOTE: mkdir for each basis set

5. For each basis set and the corresponding dir, run 

python3 00-integrals.py --basis=BASIS