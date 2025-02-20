'''
 What do we want to do?
1) Run the reconstruction software
2) Create event displays
3) Create the energy spectra plots
 '''

import os
import argparse
import numpy as np
import pandas as pd
import sys
import subprocess

def main(args):
    # define environment variables
    DATA_PATH = '/Volumes/ProtzmanSSD/data/epic/eeemcal/DESY_FEB_2025/DESY_2025/data/beam'
    OUTPUT_PATH = '/Volumes/ProtzmanSSD/data/epic/eeemcal/DESY_FEB_2025/prod'
    RUNLOG_URL = 'https://docs.google.com/spreadsheets/d/100vYwQmm6yWk3cUcB_WvoXAw8JAoOyTnIgRnm21yAfs/export?format=csv&gid=526039506'
    # WORKING_DIRECTORY = 'work'
    WORKING_DIRECTORY = '/Users/tristan/dropbox/eeemcal_desy_feb_2025'
    H2GDECODE_PATH = '/Users/tristan/epic/hgcroc/h2g_decode/build'
    ROOT_PATH = '/opt/homebrew/bin/root'

    # parse the arguments
    parser = argparse.ArgumentParser(description='Run the fast offline production')
    parser.add_argument('--run', type=int, help='Run number to process')
    parser.add_argument('--skip_decode', action='store_true', help='Skip the decoding step')

    args = parser.parse_args()
    run_number = args.run
    if run_number is None:
        print('Please provide a run number')
        return

    # create the working directory
    os.makedirs(WORKING_DIRECTORY, exist_ok=True)

    # download the runlog
    runlog = pd.read_csv(RUNLOG_URL)

    # check that the run exists in both the runlog and the data directory
    if run_number not in runlog['Run Number'].values:
        print(f'Run {run_number} not found in runlog')
        return
    h2g_file_path = os.path.join(DATA_PATH, f'Run{run_number:03}.h2g')
    if not os.path.exists(h2g_file_path):
        print(f'Run {run_number} not found in data directory')
        return
    
    # run the reconstruction software
    if not args.skip_decode:
        print(f'Running the reconstruction software on Run {run_number}')
        environment = {'DATA_PATH': DATA_PATH, 'OUTPUT_PATH': OUTPUT_PATH}
        command = [os.path.join(H2GDECODE_PATH, 'h2g_run'), str(run_number)]
        subprocess.run(command, env=environment, cwd=H2GDECODE_PATH, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('Reconstruction software finished')

    # create the event displays
    # print(f'Creating event displays for Run {run_number}')
    # command = [ROOT_PATH, '-q', '-b', '-x', '-l', f'event_display.cxx({run_number})']
    # p1 = subprocess.Popen(command, cwd=os.getcwd(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # command = [ROOT_PATH, '-q', '-b', '-x', '-l', f'event_display_tot.cxx({run_number})']
    # p2 = subprocess.Popen(command, cwd=os.getcwd(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # create the energy spectra plots
    print(f'Creating energy spectra plots for Run {run_number}')
    command = [ROOT_PATH, '-q', '-b', '-x', '-l', f'single_crystal_ADC_sum.cxx({run_number})']
    p3 = subprocess.Popen(command, cwd=os.getcwd(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #create TOT and ADC correlation plots
    print(f'Creating TOT and ADC correlation plots for Run {run_number}')
    command = [ROOT_PATH, '-q', '-b', '-x', '-l', f'adc_tot_correlation.cxx({run_number})']
    p4 = subprocess.Popen(command, cwd=os.getcwd(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # wait for the processes to finish
    # p1.wait()
    # p2.wait()
    p3.wait()
    p4.wait()

    print('Done processing, moving files...')
    os.makedirs(f'{WORKING_DIRECTORY}/run{run_number}', exist_ok=True)
    os.system(f'mv output/Run{run_number:03}*.pdf {WORKING_DIRECTORY}/run{run_number}')
    os.system(f'cp {OUTPUT_PATH}/Run{run_number:03}.root {WORKING_DIRECTORY}/run{run_number}')

    print('Fast offline production finished')
    




if __name__ == '__main__':
    main(sys.argv)