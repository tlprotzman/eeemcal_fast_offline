import pandas as pd

# EEEMCal mapping - instead of "layers", we have a single plane, where each crystal is one connector
# FPGA IP | ID
# 208     | 0
# 209     | 1
# 210     | 2
# 211     | 3
eeemcal_fpga_map = [0, 3, 3, 0, 3,
                    2, 1, 1, 1, 2,
                    2, 1, 1, 1, 3,
                    2, 2, 1, 2, 3,
                    2, 0, 0, 1, 2]

# ASIC | ID
# 0    | 0
# 1    | 1
eeemcal_asic_map = [1, 1, 1, 0, 0,
                    1, 1, 1, 1, 1,
                    1, 0, 0, 0, 0,
                    1, 0, 1, 0, 0,
                    0, 1, 1, 0, 0]

# Connector | ID
# A        | 0
# B        | 1
# C        | 2
# D        | 3
eeemcal_connector_map = [2,  0,  1,  0,  1,
                         0,  2,  0,  3,  3,
                         1,  2,  0,  3,  0,
                         2,  0,  1,  1,  2,
                         3,  1,  3,  1,  2]

eeemcal_16i_channel_a_map = [2,  6, 11, 15,  0,  4,  9, 13,
                             1,  5, 10, 14,  3,  7, 12, 16]

eeemcal_16i_channel_b_map = [20, 24, 29, 33, 18, 22, 27, 31,
                             19, 23, 28, 32, 21, 25, 30, 34]

eeemcal_16i_channel_c_map = [67, 63, 59, 55, 69, 65, 61, 57,
                             70, 66, 60, 56, 68, 64, 58, 54]

eeemcal_16i_channel_d_map = [50, 46, 40, 36, 52, 48, 42, 38,
                             51, 47, 43, 39, 49, 45, 41, 37]

eeemcal_16i_channel_map = [eeemcal_16i_channel_a_map, eeemcal_16i_channel_b_map, eeemcal_16i_channel_c_map, eeemcal_16i_channel_d_map]

eeemcal_4x4_channel_a_map = [0, 4, 9, 12]
eeemcal_4x4_channel_b_map = [19, 23, 27, 31]
eeemcal_4x4_channel_c_map = [69, 65, 61, 57]
eeemcal_4x4_channel_d_map = [52, 48, 42, 38]
eeemcal_4x4_channel_map = [eeemcal_4x4_channel_a_map, eeemcal_4x4_channel_b_map, eeemcal_4x4_channel_c_map, eeemcal_4x4_channel_d_map]

eeemcal_16p_channel_map = [6, 26, 63, 46]

sipms_per_crystal = [16, 4, 1]
crystal_ID = [5, 10, 15, 20, 25,
              4, 9, 14, 19, 24,
              3, 8, 13, 18, 23,
              2, 7, 12, 17, 22,
              1, 6, 11, 16, 21]

def make_16i_mapping():
    rows = []
    for crystal in range(25):
        for sipm in range(16):
            fpga = eeemcal_fpga_map[crystal]
            asic = eeemcal_asic_map[crystal]
            connector = eeemcal_connector_map[crystal]
            channel = eeemcal_16i_channel_map[connector][sipm] + 144 * fpga + 72 * asic
            rows.append({'FPGA': fpga, 'ASIC': asic, 'Connector': connector, 'Crystal': crystal_ID[crystal], 'SiPM': sipm, 'Channel': channel})
    df = pd.DataFrame(rows)
    df.to_csv('eeemcal_16i_mapping.csv', index=False)

def make_16p_mapping():
    rows = []
    for crystal in range(25):
        fpga = eeemcal_fpga_map[crystal]
        asic = eeemcal_asic_map[crystal]
        connector = eeemcal_connector_map[crystal]
        channel = eeemcal_16p_channel_map[connector] + 144 * fpga + 72 * asic
        rows.append({'FPGA': fpga, 'ASIC': asic, 'Connector': connector, 'Crystal': crystal_ID[crystal], 'SiPM': 1, 'Channel': channel})
    df = pd.DataFrame(rows)
    df.to_csv('eeemcal_16p_mapping.csv', index=False)

def make_4x4_mapping():
    rows = []
    for crystal in range(25):
        for sipm in range(4):
            fpga = eeemcal_fpga_map[crystal]
            asic = eeemcal_asic_map[crystal]
            connector = eeemcal_connector_map[crystal]
            channel = eeemcal_4x4_channel_map[connector][sipm] + 144 * fpga + 72 * asic
            rows.append({'FPGA': fpga, 'ASIC': asic, 'Connector': connector, 'Crystal': crystal_ID[crystal], 'SiPM': sipm, 'Channel': channel})
    df = pd.DataFrame(rows)
    df.to_csv('eeemcal_4x4_mapping.csv', index=False)

def main():
    make_16i_mapping()
    make_16p_mapping()
    make_4x4_mapping()

if __name__ == '__main__':
    main()
