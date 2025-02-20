import os
import subprocess
import pandas as pd

def main():
    OUTPUT_PATH = '/Volumes/ProtzmanSSD/data/epic/eeemcal/DESY_FEB_2025/prod'
    RUNLOG_URL = 'https://docs.google.com/spreadsheets/d/100vYwQmm6yWk3cUcB_WvoXAw8JAoOyTnIgRnm21yAfs/export?format=csv&gid=526039506'
    runlog = pd.read_csv(RUNLOG_URL)
    print(runlog.head())
    runlog = runlog[runlog['Good'] == 'GOOD']
    print(runlog.head())
    runlog['Run Number'] = runlog['Run Number'].astype(int, errors='ignore')
    runlog = runlog[runlog['Run Number'].notnull()]
    runlog = runlog[runlog['Run Number'] >= 56]
    runlog = runlog[runlog['Run Number'] <= 107]


    beam_energy_1gev_files = []

    

    for index, row in runlog.iterrows():
        # print(row['Beam Energy'])
        if row['Beam Energy'] == '4':
            print(row['Run Number'])
            file_path = os.path.join(OUTPUT_PATH, f'run{int(row['Run Number']):03}.root')
            beam_energy_1gev_files.append(file_path)

    print(beam_energy_1gev_files)

    # combine the files
    output_file = os.path.join('', 'beam_energy_1gev.root')
    command = ['/opt/homebrew/bin/hadd', output_file] + beam_energy_1gev_files
    subprocess.run(command)

if __name__ == '__main__':
    main()