import pandas as pd
import matplotlib.pyplot as plt
import datetime

# Global variables with data about line voltage and tables
mass_threads = [1]+list(range(2, 66, 2))
THREADS_AM = len(mass_threads)
LINE_VOLTAGE = 237.2
NUMBER_OF_EXP = 30
NUMBER_OF_BLOCKS = 5  # 1 original and 4 txact # Transactional module
METER_FILES = 2
RAPL_FILES = 1
LENGTH_FOR_METER_FILES = int(THREADS_AM * NUMBER_OF_EXP/METER_FILES)
LENGTH_FOR_RAPL_FILES = int(THREADS_AM * NUMBER_OF_EXP/RAPL_FILES)

# DataMeters and DataMetersTime
dm = pd.read_csv('dataMeters.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt = pd.read_csv('dataMetersTime.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])

# DataMeters2 and DataMetersTime2
dm2 = pd.read_csv('dataMeters2.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt2 = pd.read_csv('dataMetersTime2.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])

# DataRAPL
dr = pd.read_csv('dataRAPL.csv', sep=',', usecols=['version', 'n-workers', 'n-secondary-workers', 'n-reservations', 'n-relations', 'n-queries', 'password-work-factor', 'Energy Cores', 'Energy Ram', 'Energy Gpu', 'Energy Pkg', 'Time (perf) (s)', 'Time (exec) (s)', 'Ratio time (%)', 'Power_Cores', 'Power_Ram', 'Power_Gpu', 'Power_Pkg'])

# DataDELTAS
dD = pd.read_csv('dataDELTAS.csv', sep=',', usecols=['version', 'n-workers', 'n-secondary-workers', 'n-reservations', 'n-relations', 'n-queries', 'password-work-factor', 'Energy Delta'])

# General min-max
print(max(dD["Energy Delta"]))  # 295.23071600000026
print(min(dD["Energy Delta"]))  # 0.0270519999999407

# Original min-max
print(max(dD["Energy Delta"][:990]))  # 193.29
print(min(dD["Energy Delta"][:990]))  # 0.1427760000002251

def get_zip_startend_time(source, ind):
    """ Returns the time frames from the source table

    @:param source: The used table
    @:type source: pandas type
    @:param ind: The index of the part with time
    @:type ind: int
    @:returns: A zipped set with start and end time
    @:rtype: set
    """
    starts = [i.split(" ")[ind] for i in source["Start Time"]]
    ends = starts[1:] + [(source["Stop Time"][-1:].values[0]).split(" ")[ind]]
    return zip(starts, ends)


def get_data_for_15(a, b, source_dm):
    """ Calculates, with the provided formulas, averages for the 1...64 threads
    in first 15 blocks and then returns the list with them

    formula: Energy[n] =  C9[n] x U x C5[n]

    Where:
    C9 is the Average
    U is the line voltage, which was equal to 237.2 Volts
    C5[n] is the duration

    after getting the list, divide its sum by 15 (blocks)

    @:param a: The list with time frames to be used
    @:type a: set
    @:param b: The list with time frames to be checked
    @:type b: int
    @:param source_dm: A DataMeters file to be used
    @:type source_dm: pandas type
    @:returns: a list with averages for the 1...64 threads in first fifteen blocks
    @:rtype: list
    """
    id_dict, energy_ress, mega_b, cc = 1, {i: [] for i in mass_threads}, b, 1
    ## id_dict, energy_ress, mega_b, cc, add_l = 1, {i: [] for i in mass_threads}, b, 1,[] #  wrong solution becase average frames
    for (start_a, end_a) in a:
        # print(">>> ",cc," <<< ") if id_dict ==1 else None
        # print(">>",id_dict, " <<")
        # print("    ", start_a," --- ", end_a)
        for (start_b, end_b) in b:
            if start_a <= start_b and end_b <= end_a:
                energy_ress[id_dict] +=[float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',','.')) * LINE_VOLTAGE * float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',','.'))]
                ## add_l += [float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',','.')) * LINE_VOLTAGE * float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',','.'))]
                # print("    " * 2, start_b, " -- ", end_b, " ~ ", float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',', '.')), " * ", LINE_VOLTAGE, " * ", float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',', '.')))
            if start_b > start_a and end_b > end_a:
                ## res = sum(add_l)/len(add_l) if len(add_l)!=0 else 0
                ## energy_ress[id_dict]+=[res]
                ## add_l=[]
                b = b[b.index((start_b, end_b))-1:]
                break
        id_dict += 1 if id_dict == 1 else 2
        if id_dict == 66: id_dict = 1; cc += 1
    # return [sum(energy_ress[i])/len(energy_ress[i]) for i in mass_threads]
    return [sum(energy_ress[i])/15 for i in mass_threads]


def get_time_in_time(dm_tu, dmt_tu):
    """ Returns the set with time frame in time format

    @:param dm_tu: The used DataMeters table
    @:type dm_tu: pandas type
    @:param dmt_tu: The used DataMetersTime table
    @:type dmt_tu: pandas type
    @:returns: A set with start and end time in time format
    @:rtype: set
    """
    a_dt = [(datetime.datetime.strptime(start, '%H:%M:%S'), datetime.datetime.strptime(end, '%H:%M:%S')) for start, end in get_zip_startend_time(dmt_tu,4)]
    b_dt = [(datetime.datetime.strptime(start, '%H:%M:%S,%f'),datetime.datetime.strptime(end, '%H:%M:%S,%f')) for start, end in get_zip_startend_time(dm_tu,1)]
    return (a_dt, b_dt)


def get_res_data_rapl(dr_tu):
    """ Returns a list with all Energy consumption related columns row sums

    @:param dr_tu: The used DataRAPL table
    @:type dr_tu: pandas type
    @:returns: A list with sums
    @:rtype: list
    """
    # return [i + j + k + t for (i, j, k, t) in zip(dr_tu['Energy Cores'], dr_tu['Energy Ram'], dr_tu['Energy Gpu'], dr_tu['Energy Pkg'])]
    # return [i + j + k  for (i, j, k) in zip(dr_tu['Energy Cores'], dr_tu['Energy Ram'], dr_tu['Energy Gpu'])] # checking without Energy Pkg
    return [j + k + t for (j, k, t) in zip( dr_tu['Energy Ram'], dr_tu['Energy Gpu'], dr_tu['Energy Pkg'])] # checking without cores


def calc_avar(source):
    """ Returns a list with averages for the 1...64 threads for RAPL

    @:param source: The list with data to be used
    @:type source: list
    @:returns: A list with sums for each thread
    @:rtype: list
    """
    mass = []
    for i in range(THREADS_AM): mass.append(sum(list(source[i::THREADS_AM])) / ((len(source) / THREADS_AM)))
    return mass

# CALLING FUNCTIONS:

# getting time frames in time format for DatMeters and DatMetersTime tables
(a_dt_1, b_dt_1) = get_time_in_time(dm, dmt)
(a_dt_2, b_dt_2) = get_time_in_time(dm2, dmt2)

# getting the averages for the 15 blocks for two files and
# saving them into a dictionary with key saying if it is original(1) or txact(>1)
mass_with_ress_Meter = {(i+1): (get_data_for_15(a_dt_1[LENGTH_FOR_METER_FILES*i:LENGTH_FOR_METER_FILES*(i+1)], b_dt_1, dm), get_data_for_15(a_dt_2[LENGTH_FOR_METER_FILES * i:LENGTH_FOR_METER_FILES *(i+1)], b_dt_2, dm2)) for i in range(5)}

# because we need average for 30 blocks just summing respectful results and divide by 2 (15+15 = 30 blocks)
dict_with_Meter = {}
for j in range(1, NUMBER_OF_BLOCKS+1):
    dict_with_Meter[j] = [(mass_with_ress_Meter[j][0][i] + mass_with_ress_Meter[j][1][i])/2 for i in range(THREADS_AM)]

# getting data from RAPL and calculating averages for all 30 blocks (5)
mass_with_ress_RAPL = get_res_data_rapl(dr)

dict_with_RAPL = {}
for j in range(1,NUMBER_OF_BLOCKS+1):
    dict_with_RAPL[j] = calc_avar(mass_with_ress_RAPL[LENGTH_FOR_RAPL_FILES*(j-1):LENGTH_FOR_RAPL_FILES*j])

# Creating plots with all collected data
figure, axis = plt.subplots(2)

# METERS PLOT (TOP)
axis[0].plot(mass_threads, dict_with_Meter[1], "b", label="original")
axis[0].plot(mass_threads, dict_with_Meter[2], "r", label="txact s1")
axis[0].plot(mass_threads, dict_with_Meter[3], "g", label="txact s2")
axis[0].plot(mass_threads, dict_with_Meter[4], "c", label="txact s8")
axis[0].plot(mass_threads, dict_with_Meter[5], "m", label="txact s64")
axis[0].set_title("Meters data (origin blue) energy consumption 30")
# axis[0].set_xlabel('1-64 threads')
axis[0].set_ylabel('energy consumption (J)')
axis[0].legend()

# RAPL PLOT (BOTTOM)
axis[1].plot(mass_threads, dict_with_RAPL[1], "b", label="original")
axis[1].plot(mass_threads, dict_with_RAPL[2], "r", label="txact s1")
axis[1].plot(mass_threads, dict_with_RAPL[3], "g", label="txact s2")
axis[1].plot(mass_threads, dict_with_RAPL[4], "c", label="txact s8")
axis[1].plot(mass_threads, dict_with_RAPL[5], "m", label="txact s64")
# axis[1].set_title("RAPL data (origin magenta) energy consumtion 30")
axis[1].set_xlabel('RAPL data (origin magenta) energy consumption 30')
axis[1].set_ylabel('energy consumption (J)')

axis[1].legend()
plt.show()
