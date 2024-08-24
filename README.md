# pPb "UPC" skim code using HTCondor

Code to produce jets, tracks and event plane skims from the CMS HiForest. Setup here: https://github.com/denerslemos/JetSmallSystems/tree/pPb/HiForest

## Intructions

Setup CMSSW (just for root versioning)
```
export SCRAM_ARCH=el9_amd64_gcc12
cmsrel CMSSW_13_0_5
cd CMSSW_13_0_5/src
cmsenv
```
Inside of the src folder, download the code using
```
git clone https://github.com/denerslemos/pPbSkimsUPC.git
cd pPbSkimsUPC
mkdir cond
```
Before compile the code you must check the [sub_skim.sh](https://github.com/denerslemos/pPbSkimsUPC/blob/main/sub_skim.sh) lines 4 (CMSSW/src) and 6 (.../pPb2016skims) and replace by your own folders.

Once this steps are done you can compile the code with
```
g++ -O2 pPbSkimsUPC.C `root-config --libs` `root-config --cflags` -o pPbSkimsUPC
```
This will create the executable: ```pPbSkimsUPC``` 

After that you will need your VOMS certificate, do it using
```
voms-proxy-init -rfc -voms cms --out voms_proxy.txt --hours 200
```
that creates a certificate file valid for 200 hours: ```voms_proxy.txt```

Now you can submit the condor jobs using the python script, [```HTCondor_submit.py```](https://github.com/denerslemos/pPbSkimsUPC/blob/main/HTCondor_submit.py):

```
python3 HTCondor_submit.py -i input_text_file -o output_name_file -m X -n Y -s Z -x Q -z W
```

- input_text_file: is the text file (use it without the .txt extension) with inputs and can be found in the folders [MC_PYTHIA_SAMPLES](https://github.com/denerslemos/pPbSkimsUPC/tree/main/MC_PYTHIA_SAMPLES) or [DATA_SAMPLES](https://github.com/denerslemos/pPbSkimsUPC/tree/main/DATA_SAMPLES) each .root input will be a job

- output_name_file: output file name (use it without the .root extension), it will automatically include a counter for each input. You can use paths to save on EOS.

- X: 0 for data and 1 for MC

- Y: 0 for no multiplicity cut (mostly MC or jet samples), 1 for MB [0,185], 2 for HM185 [185,250] and 3 for HM250 [250,inf]

- Z: name for the submission files, I have used HTcondor_sub_ + some information from the sample, pthat, MB, ... + pgoing or Pbgoing.

- Q: data side MB triggers due  to different names : 0 for p -> + eta and 1 for p -> - eta. For MC use always 0.

- W: Files with ZDC information, it depends on the MB PD.


It will automatically include a counter for each input
