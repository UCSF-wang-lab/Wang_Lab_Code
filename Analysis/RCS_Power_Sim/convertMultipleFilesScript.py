import glob
import subprocess

data_path = "/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS09/v6_2024-02-28_aDBS-001/Data/Analysis Data/aDBS/ConvertedData"
save_path = data_path

left_fftsize = "256"
left_bitshift = "1"
left_files = glob.glob(data_path+"/*LEFT.csv")

for f in left_files:
    subprocess.run(["python3","calculateRCSPower.py","-data",f,"-settings",f.replace("LEFT","LEFT_SETTINGS"),"-nfft",left_fftsize,"-fft_bitshift",left_bitshift,"-save_path",save_path])


right_fftsize = "256"
right_bitshift = "2"
right_files = glob.glob(data_path+"/*RIGHT.csv")

for f in right_files:
    subprocess.run(["python3","calculateRCSPower.py","-data",f,"-settings",f.replace("RIGHT","RIGHT_SETTINGS"),"-nfft",right_fftsize,"-fft_bitshift",right_bitshift,"-save_path",save_path])
