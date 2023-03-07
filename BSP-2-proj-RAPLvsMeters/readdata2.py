import pandas as pd
import matplotlib.pyplot as plt
import datetime

# Global variables with data from tables
mass_threads= [1]+list(range(2, 66, 2))

# DataMeters and DataMetersTime
dm = pd.read_csv('dataMeters.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt = pd.read_csv('dataMetersTime.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])

# DataMeters2 and DataMetersTime2
dm2 = pd.read_csv('dataMeters2.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt2 = pd.read_csv('dataMetersTime2.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])

# DataRAPL
dr = pd.read_csv('dataRAPL.csv', sep=',', usecols=['version', 'n-workers', 'n-secondary-workers', 'n-reservations', 'n-relations', 'n-queries', 'password-work-factor', 'Energy Cores', 'Energy Ram', 'Energy Gpu', 'Energy Pkg', 'Time (perf) (s)', 'Time (exec) (s)', 'Ratio time (%)', 'Power_Cores', 'Power_Ram', 'Power_Gpu', 'Power_Pkg'])


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
    for (start_a, end_a) in a:
        # print(">>> ",cc," <<< ") if id_dict ==1 else None
        # print(">>",id_dict, " <<")
        # print("    ", start_a," --- ", end_a)
        for (start_b, end_b) in b:
            if start_a <= start_b and end_b <= end_a:
                energy_ress[id_dict] +=[float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',','.')) * 237.2 * float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',','.'))]
                # print(float(str(source_dm["Duration"][b.index((start_b, end_b))]).split(":")[2].replace(',','.')))
                # print("    " * 2, start_b, " -- ", end_b, " ~ ", float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',', '.')), " * ", 237.2, " * ", float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',', '.')))
            if start_b > start_a and end_b > end_a:
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
    return [i + j + k + t for (i, j, k, t) in zip(dr_tu['Energy Cores'], dr_tu['Energy Ram'], dr_tu['Energy Gpu'], dr_tu['Energy Pkg'])]


def calc_avar(mass, source):
    """ Adds, to the passed to the function list, averages for the 1...64 threads

    @:param mass: The list to be changed
    @:type mass: list
    @:param source: The list with data to be used
    @:type source: list
    """
    for j in range(33): mass.append(sum(list(source[j::33])) / ((len(source) / 33)))


# CALLING FUNCTIONS:

# getting time frames in time format for DatMeters and DatMetersTime tables
(a_dt_1,b_dt_1) = get_time_in_time(dm, dmt)
(a_dt_2,b_dt_2) = get_time_in_time(dm2, dmt2)

# getting the averages for the 15 blocks for two files and
# saving them into a dictionary with key saying if it is original(1) or txact(>1)
mega_mass_with_mega_data = {(i+1): (get_data_for_15(a_dt_1[495*i:495+495*i], b_dt_1, dm), get_data_for_15(a_dt_2[495 * i:495 + 495 * i], b_dt_2, dm2)) for i in range(5)}

# because we need average for 30 blocks just summing respectful results and divide by 2 (15+15 = 30 blocks)
newMMWMD = {}
for j in range(1, 6):
    newMMWMD[j] = [(mega_mass_with_mega_data[j][0][i] + mega_mass_with_mega_data[j][1][i])/2 for i in range(33)]

# getting data from RAPL and calculating averages for all 30 blocks (5)
mass_with_ress = get_res_data_rapl(dr)

# print(len(mass_with_ress))  #4950
mass_data_rapl_or =[]  # original
mass_data_rapl_tx1 =[]  # txact
mass_data_rapl_tx2 =[]  # txact
mass_data_rapl_tx3 =[]  # txact
mass_data_rapl_tx4 =[]  # txact
calc_avar(mass_data_rapl_or,mass_with_ress[0:990])
calc_avar(mass_data_rapl_tx1,mass_with_ress[990:1980])
calc_avar(mass_data_rapl_tx2,mass_with_ress[1980:2970])
calc_avar(mass_data_rapl_tx3,mass_with_ress[2970:3960])
calc_avar(mass_data_rapl_tx4,mass_with_ress[3960:4950])

# Creating plots with all collected data
figure, axis = plt.subplots(2)

# METERS PLOT (TOP)
axis[0].plot(mass_threads, newMMWMD[1], "b", label= "origin")
axis[0].plot(mass_threads, newMMWMD[2], "r", label= "txact")
axis[0].plot(mass_threads, newMMWMD[3], "r", label= "txact")
axis[0].plot(mass_threads, newMMWMD[4], "r", label= "txact")
axis[0].plot(mass_threads, newMMWMD[5], "r", label= "txact")
axis[0].set_title("Meters data (origin blue) energy consumtion 30")
# axis[0].set_xlabel('1-64 threads')
axis[0].set_ylabel('energy consumption (J)')
axis[0].legend()

# RAPL PLOT (BOTTOM)
axis[1].plot(mass_threads, mass_data_rapl_or, "m", label= "origin")
axis[1].plot(mass_threads, mass_data_rapl_tx1, "c", label= "txact")
axis[1].plot(mass_threads, mass_data_rapl_tx2, "c", label= "txact")
axis[1].plot(mass_threads, mass_data_rapl_tx3, "c", label= "txact")
axis[1].plot(mass_threads, mass_data_rapl_tx4, "c", label= "txact")
# axis[1].set_title("RAPL data (origin magenta) energy consumtion 30")
axis[1].set_xlabel('RAPL data (origin magenta) energy consumtion 30')
axis[1].set_ylabel('energy consumption (J)')

axis[1].legend()
plt.show()
