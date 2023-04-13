import pandas as pd
import matplotlib.pyplot as plt
import datetime
import glob
import os
from colorama import Fore

# Global variables with data about line voltage and tables
n_secondary_workers = [0,1,2,8,64]
mass_threads = [1]+list(range(2, 66, 2))
THREADS_AM = len(mass_threads)
LINE_VOLTAGE = 237.2
NUMBER_OF_EXP = 30
TRANS_MODEL = 5  # 1 original and 4 txact
METER_FILES = 2
RAPL_FILES = 1
LENGTH_FOR_METER_FILES = int(THREADS_AM * NUMBER_OF_EXP/METER_FILES)
LENGTH_FOR_RAPL_FILES = int(THREADS_AM * NUMBER_OF_EXP/RAPL_FILES)


def get_columns_for_file(filename, sep):
    """ Gets the column names for a certain file

    @:param filename: name of the file
    @:type filename: string
    @:param sep: separator between words
    @:type sep: string
    @:return mass: list of column names
    @:rtype: list
    """
    with open(filename, 'r', encoding='utf-8') as file:
        mass = file.readline().split(sep)
        mass[-1]= mass[-1][:-1]
        mass = list(filter(lambda x: x!="",mass))
        return mass


def get_sep(filename):
    """ Gets the separator for a certain file

    @:param filename: name of the file
    @:type filename: string
    @:return i: separator
    @:rtype: string
    """
    with open(filename, 'r', encoding='utf-8') as file:
        for i in file.readline():
            if not i.isalpha() and i != " ":
                return i


def get_lists_of_combo_Meters():
    """ Creates a list of tuples with file and its separator (for Time and Meter files)

    @:return megalist: the list of tuples with file and its separator
    @:rtype: list
    """
    listfiles,dms_l,dmts_l=[],[],[]
    for i in glob.glob(os.path.abspath('*.csv')): listfiles += [i.split("\\")[-1]]
    listMeterfiles = list(filter(lambda x: x.startswith("dataMeters"), listfiles))
    for i in listMeterfiles:
        if i.startswith("dataMetersTime"): dmts_l.append(i)
        else: dms_l.append(i)
    megalist = []
    for (t,m) in zip(dmts_l,dms_l):  megalist+=((t,get_sep(t)),(m,get_sep(m)))
    return megalist


def get_pandas_files_Meter():
    """ Creates a list of files in pandas type

    @:return pandas_f: list of files in pandas type
    @:rtype: list
    """
    files_with_seps = get_lists_of_combo_Meters()
    pandas_f=[]
    for (f,s) in files_with_seps:
        usecols_l = get_columns_for_file(f, s)
        pandas_f.append(pd.read_csv(f, sep=s, usecols=usecols_l))
    return pandas_f

list_with_pandas_files_meter = get_pandas_files_Meter()

# DataMeters and DataMetersTime
dmt_tup = tuple()
dm_tup=tuple()
for i in range(int(len(list_with_pandas_files_meter)/2)):
    dmt_tup = dmt_tup+ (list_with_pandas_files_meter[2*i],)
    dm_tup = dm_tup+ (list_with_pandas_files_meter[2*i+1],)

sep_for_RMD = get_sep('dataRAPL.csv')
# DataRAPL
dR = pd.read_csv('dataRAPL.csv', sep=sep_for_RMD, usecols=get_columns_for_file('dataRAPL.csv',sep_for_RMD))
dM = pd.read_csv('dataMETER.csv', sep=sep_for_RMD,usecols=get_columns_for_file('dataMETER.csv',sep_for_RMD) )
dD = pd.read_csv('dataDELTAS.csv', sep=sep_for_RMD,usecols=get_columns_for_file('dataDELTAS.csv',sep_for_RMD))

dM_v2 = pd.read_csv('dataMETER_ver2.csv', sep=sep_for_RMD,usecols=get_columns_for_file('dataMETER_ver2.csv',sep_for_RMD) )
dD_v2 = pd.read_csv('dataDELTAS_ver2.csv', sep=sep_for_RMD,usecols=get_columns_for_file('dataDELTAS_ver2.csv',sep_for_RMD))


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


def get_data_for_15_ver_2(a, b, source_dm ,rewritefile = False, mid_val_th=tuple(),approach=False):
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
    @:param rewritefile: The flag which indicates if I am trying to get results for each thread amount
    @:type rewritefile: boolean
    @:param mid_val_th: The tuple which contains energy consumption for each row
    @:type mid_val_th: tuple
    @:returns: a list with averages for the 1...64 threads in first fifteen blocks if rewritefile is false, else mid_val_th with energy consumptions
    @:rtype: list / tuple
    """
    id_dict, energy_ress, mega_b, cc,mid_val_th_calc = 1, {i: [] for i in mass_threads}, b, 1,[]

    ## id_dict, energy_ress, mega_b, cc, add_l = 1, {i: [] for i in mass_threads}, b, 1,[] #  wrong solution becase average frames
    for (start_a, end_a) in a:
        print(">>> ",cc," <<< ") if id_dict ==1 else None
        print(">>",id_dict, " <<")
        print("    ", start_a," --- ", end_a)
        for (start_b, end_b) in b:
            if (approach and start_a <= start_b and start_b < end_a) or start_a <= start_b and end_b <= end_a:
                if (rewritefile): mid_val_th_calc += [float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',','.')) * LINE_VOLTAGE * float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',','.'))]
                else: energy_ress[id_dict] +=[float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',','.')) * LINE_VOLTAGE * float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',','.'))]
                ## add_l += [float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',','.')) * LINE_VOLTAGE * float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',','.'))]
                print("    " * 2, start_b, " -- ", end_b, " ~ ", float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',', '.')), " * ", LINE_VOLTAGE, " * ", float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',', '.')))
            if (approach and start_b >= end_a)   or start_b > start_a and end_b > end_a:
                ## res = sum(add_l)/len(add_l) if len(add_l)!=0 else 0
                ## energy_ress[id_dict]+=[res]
                ## add_l=[]
                ## b = b[b.index((start_b, end_b))-1:] check if -1 is needed
                b = b[b.index((start_b, end_b)):] ## no difference
                break
        if (rewritefile):
            mid_val_th = mid_val_th + (sum(mid_val_th_calc),)
            mid_val_th_calc=[]
            continue
        id_dict += 1 if id_dict == 1 else 2
        if id_dict == 66: id_dict = 1; cc += 1
    # return [sum(energy_ress[i])/len(energy_ress[i]) for i in mass_threads]
    if (rewritefile): return mid_val_th
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
    ## return [i + j + k + t for (i, j, k, t) in zip(dr_tu['Energy Cores'], dr_tu['Energy Ram'], dr_tu['Energy Gpu'], dr_tu['Energy Pkg'])]
    ## return [i + j + k  for (i, j, k) in zip(dr_tu['Energy Cores'], dr_tu['Energy Ram'], dr_tu['Energy Gpu'])]  # checking without Energy Pkg
    return [j + k + t for (j, k, t) in zip(dr_tu['Energy Ram'], dr_tu['Energy Gpu'], dr_tu['Energy Pkg'])]  # checking without cores


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


def rewriteFileMETERandDeltas(filename1, filename2, boolV2):
    """ Rewrites data in files dataMETER and dataRAPL """
    # getting time frames in time format for DatMeters and DatMetersTime tables
    (a_dt_1, b_dt_1) = get_time_in_time(dm_tup[0], dmt_tup[0])
    (a_dt_2, b_dt_2) = get_time_in_time(dm_tup[1], dmt_tup[1])
    # getting data from RAPL file to get deltas
    dr_data = get_res_data_rapl(dR)

    with open(filename1, 'w', encoding='utf-8') as file:
        with open(filename2, 'w', encoding='utf-8') as file2:
            file.write("version,n-workers,n-secondary-workers,n-reservations,n-relations,n-queries,password-work-factor,Energy Meter\n")
            file2.write("version,n-workers,n-secondary-workers,n-reservations,n-relations,n-queries,password-work-factor,Energy Delta\n")

    newtuple=tuple()
    for i in range(TRANS_MODEL):
        newtuple = get_data_for_15_ver_2(a_dt_1[LENGTH_FOR_METER_FILES*i:LENGTH_FOR_METER_FILES*(i+1)], b_dt_1,dm_tup[0], True, newtuple,boolV2)
        newtuple = get_data_for_15_ver_2(a_dt_2[LENGTH_FOR_METER_FILES*i:LENGTH_FOR_METER_FILES*(i+1)], b_dt_2, dm_tup[1], True, newtuple,boolV2)

    mt = 1
    with open(filename1, 'a', encoding='utf-8') as file:
        with open(filename2, 'a', encoding='utf-8') as file2:
            for ind_for_mi,nsw in enumerate(n_secondary_workers):
                for mi in range(ind_for_mi*LENGTH_FOR_METER_FILES*2,(ind_for_mi+1)* LENGTH_FOR_METER_FILES*2):
                    energyM = newtuple[mi]
                    file.write(f"{'original' if nsw == 0 else 'txact'},{mt},{nsw},1000,50,10,5,{energyM}\n")
                    file2.write(f"{'original' if nsw == 0 else 'txact'},{mt},{nsw},1000,50,10,5,{abs(energyM - dr_data[mi])  }\n")
                    mt += 1 if mt == 1 else 2
                    mt =1 if mt not in mass_threads else mt


def get_distribution(dD_frame):
    """ Creates a dictionary with amount of energies for each frame (10>,10-50,...,200<)

    @:param dD_frame: the frame to be checked
    @:type dD_frame: pandas
    @:return dict_with_disbs: the dictionary with amount of energies for each frame
    @:rtype: dictionary
    """
    dict_with_disbs={"10 >":[],"10-50":[],"50-100":[],"100-200":[],"200 <":[]}
    for i in dD_frame:
        if i<10: dict_with_disbs["10 >"]+=[i]
        if i>10 and i <50: dict_with_disbs["10-50"]+=[i]
        if i>=50 and i <100: dict_with_disbs["50-100"]+=[i]
        if i>=100 and i <200: dict_with_disbs["100-200"]+=[i]
        if i>=200: dict_with_disbs["200 <"]+=[i]
    return dict_with_disbs


# RMDe = [get_res_data_rapl(dR),dM["Energy Meter"],dD["Energy Delta"]]
RMDe_2 = [get_res_data_rapl(dR),dM["Energy Meter"],dD["Energy Delta"],get_res_data_rapl(dR),dM_v2["Energy Meter"],dD_v2["Energy Delta"]]

def print_for_each_in_RMD(frame):
    """ Prints min, max and distribution for a certain frame

        @:param frame: the list with first and last indices
        @:type frame: list
        """
    datamass = [i[frame[0]:frame[1]] for i in RMDe_2]
    print("\n"," "*19,"RAPL"," "*16,"METER"," "*15,"DELTA", end="")
    print("        ", " " * 19, "RAPL", " " * 16, "METER_2", " " * 15, "DELTA_2:")

    print("max:     ", end="")
    for c,i in enumerate(datamass):
        if ((c+1)%4)==0:print('{:>14}'.format(" "), end="")
        print('{:>22}'.format(max(i)), end="")

    print("\nmin:     ", end="")
    for c,i in enumerate(datamass):
        if ((c+1)%4)==0: print('{:>14}'.format(" "), end="")
        print('{:>22}'.format(min(i)), end="")
    print()

    dict_with_disbs = [get_distribution(i) for i in datamass]
    for i in ["10 >", "10-50", "50-100", "100-200", "200 <"]:
        print('\n{:<9}'.format(f"{i}: "),end="")
        for c,j in enumerate(dict_with_disbs):
            if ((c+1)%4)==0: print('{:>14}'.format(" "), end="")
            print('{:>22}'.format( len(j[i])),end="")

    means = [sum(i)/len(i) for i in datamass]
    print("\n\nMean:    ",end="")
    for c,i in enumerate(means):
        if ((c+1)%4)==0: print('{:>14}'.format(" "), end="")
        print('{:>22}'.format(i), end="")

    deviations = [sum([pow((i-means[j]),2) for i in datamass[j]])/(LENGTH_FOR_METER_FILES*2) for j in range(len(datamass)) ]
    print("\nSt. Dev: ", end="")
    for c,i in enumerate(deviations):
        if ((c+1)%4)==0: print('{:>14}'.format(" "), end="")
        print('{:>22}'.format(pow(i,0.5)), end="")


def get_everything_about_RMD():
    """ Prints min, max and distribution for RAPL, Meter and Deltas for original and txacts """
    colors,space = [Fore.BLUE, Fore.RED, Fore.GREEN, Fore.CYAN, Fore.MAGENTA,Fore.YELLOW],138
    print(colors[0] + " --- ORIGINAL: ","-"*space)
    print_for_each_in_RMD([0,LENGTH_FOR_METER_FILES*2])
    for n in range(len(n_secondary_workers)-1):
        print(colors[n+1] + f"\n\n --- TXACT-{n_secondary_workers[n+1]}:  ","-"*space)
        print_for_each_in_RMD([LENGTH_FOR_METER_FILES*(2)*(n+1), LENGTH_FOR_METER_FILES * (2)*(n+2)])

    print(colors[-1]+ "\n\n --- GENERAL:   ","-"*space)
    print_for_each_in_RMD([0,LENGTH_FOR_METER_FILES*2*5])


def show_plot(flag_for_ver_2):
    """ Shows graphical representation of data from RAPL and METER"""
    # getting the averages for the 15 blocks for two files and
    # saving them into a dictionary with key saying if it is original(1) or txact(>1)
    # mass_with_ress_Meter = {(i+1): (get_data_for_15(a_dt_1[LENGTH_FOR_METER_FILES*i:LENGTH_FOR_METER_FILES*(i+1)], b_dt_1, dm), get_data_for_15(a_dt_2[LENGTH_FOR_METER_FILES * i:LENGTH_FOR_METER_FILES *(i+1)], b_dt_2, dm2)) for i in range(5)}
    #
    # # because we need average for 30 blocks just summing respectful results and divide by 2 (15+15 = 30 blocks)
    # dict_with_Meter = {}
    # for j in range(1, TRANS_MODEL+1):
    #     dict_with_Meter[j] = [(mass_with_ress_Meter[j][0][i] + mass_with_ress_Meter[j][1][i])/2 for i in range(THREADS_AM)]

    ## example with data after calculations (problems with average)
    ## mass_with_ress_METER = [i for i in dM["Energy Meter"]]
    ##
    ## dict_with_Meter = {}
    ## for j in range(1,TRANS_MODEL+1):
    ##     dict_with_Meter[j] = calc_avar(mass_with_ress_METER[LENGTH_FOR_RAPL_FILES*(j-1):LENGTH_FOR_RAPL_FILES*j])
    ##
    ## Result Nothing changed

    # getting data from RAPL and Meter and calculating averages for all 30 blocks (5)

    if flag_for_ver_2: mass_with_ress_RAPL,mass_with_ress_METER = get_res_data_rapl(dR),[i for i in dM_v2["Energy Meter"]]
    else:mass_with_ress_RAPL,mass_with_ress_METER = get_res_data_rapl(dR),[i for i in dM["Energy Meter"]]

    dict_with_Meter,dict_with_RAPL = {},{}
    for j in range(1,TRANS_MODEL+1):
        dict_with_RAPL[j] = calc_avar(mass_with_ress_RAPL[LENGTH_FOR_RAPL_FILES*(j-1):LENGTH_FOR_RAPL_FILES*j])
        dict_with_Meter[j] = calc_avar(mass_with_ress_METER[LENGTH_FOR_RAPL_FILES * (j - 1):LENGTH_FOR_RAPL_FILES * j])

    # Creating plots with all collected data
    figure, axis = plt.subplots(2)

    # METERS PLOT (TOP)
    # RAPL PLOT (BOTTOM)
    # colors
    colors= ['b','r','g','c','m','indigo','k','y','darkorange']

    axis[0].plot(mass_threads, dict_with_Meter[1], colors[0], label="original")
    axis[1].plot(mass_threads, dict_with_RAPL[1], colors[0], label="original")

    for i in range(1,TRANS_MODEL):
        axis[0].plot(mass_threads, dict_with_Meter[i+1], colors[i], label=f"txact s{n_secondary_workers[i]}")
        axis[1].plot(mass_threads, dict_with_RAPL[i+1], colors[i], label=f"txact s{n_secondary_workers[i]}")

    axis[0].set_title("Meters data (origin blue) energy consumption 30")
    # axis[0].set_xlabel('1-64 threads')
    axis[0].set_ylabel('energy consumption (J)')
    axis[0].legend()

    # axis[1].set_title("RAPL data (origin magenta) energy consumtion 30")
    axis[1].set_xlabel('RAPL data (origin magenta) energy consumption 30')
    axis[1].set_ylabel('energy consumption (J)')

    axis[1].legend()
    plt.show()


def get_plot_orig_time_rapl_avar():
    """ Shows speed-up for each number of threads in RAPL using averages """
    dict_with_RAPL,txacts = {},{}
    for j in range(1, TRANS_MODEL + 1): dict_with_RAPL[j] = calc_avar(list(dR["Time (exec) (s)"][LENGTH_FOR_RAPL_FILES * (j - 1): LENGTH_FOR_RAPL_FILES * j]))
    orig_mass=[i/dict_with_RAPL[1][0] for i in dict_with_RAPL[1]]
    for t in range(1, TRANS_MODEL): txacts[t] = [(j/i)*c  for (i,j,c) in zip(dict_with_RAPL[t+1],dict_with_RAPL[1], orig_mass)]

    colors = ['b', 'r', 'g', 'c', 'm', 'indigo', 'k', 'y', 'darkorange']
    plt.plot(mass_threads,orig_mass,colors[0],label=f"original")
    for i in range(1, len(txacts)+1):
        plt.plot(mass_threads,txacts[i],colors[i],label=f"txact-{n_secondary_workers[i]}")
    plt.xlabel('number of threads')
    plt.ylabel('Speed-up')
    plt.title("RAPL speed-up check")
    plt.legend()
    plt.show()

def get_plot_orig_time_rapl_median():
    """ Shows speed-up for each number of threads in RAPL using medians """
    dict_with_RAPL,txacts = {i:[] for i in range(1,TRANS_MODEL+1)},{}
    for j in range(1, TRANS_MODEL + 1):
        mass = list(dR["Time (exec) (s)"][LENGTH_FOR_RAPL_FILES * (j - 1): LENGTH_FOR_RAPL_FILES * j])
        for p in range(THREADS_AM):
            f = mass[p::33]
            f.sort()
            dict_with_RAPL[j] += [f[int((len(f)-1)/2)]]

    orig_mass=[i/dict_with_RAPL[1][0] for i in dict_with_RAPL[1]]
    for t in range(1, TRANS_MODEL): txacts[t] = [(j/i)*c  for (i,j,c) in zip(dict_with_RAPL[t+1],dict_with_RAPL[1], orig_mass)]

    colors = ['b', 'r', 'g', 'c', 'm', 'indigo', 'k', 'y', 'darkorange']
    plt.plot(mass_threads,orig_mass,colors[0],label=f"original")
    for i in range(1, len(txacts)+1):
        plt.plot(mass_threads,txacts[i],colors[i],label=f"txact-{n_secondary_workers[i]}")
    plt.xlabel('number of threads')
    plt.ylabel('Speed-up')
    plt.title("RAPL speed-up check")
    plt.legend()
    plt.show()


def check_data_before_meter_start(dmt,dM_data):
    """ Calculates average energy consumption before the start of the experiment

            @:param dmt: table with all Meters execution description
            @:type dmt: pandas type
            @:param dM_data: table with all Meters execution results
            @:type dM_data: pandas type
            @returns average energy consumption obtained before the start of the experiment
        """
    timelim1,list_with_st_time,energy_consumed = datetime.datetime.strptime( dmt["Start Time"][0].split()[3], '%H:%M:%S'),[],[]
    for i in dM_data.loc:
        if datetime.datetime.strptime(i["Stop Time"].split()[1], '%H:%M:%S,%f') > timelim1:
            break
        energy_consumed += [float(str(i["Average"]).replace(',', '.')) * LINE_VOLTAGE * float(str(i["Duration"]).split(":")[2].replace(',', '.'))]
    return sum(energy_consumed)/len(energy_consumed)

def check_data_after_meter_start(dmt,dM_data):
    """ Calculates average energy consumption after the start of the experiment

        @:param dmt: table with all Meters execution description
        @:type dmt: pandas type
        @:param dM_data: table with all Meters execution results
        @:type dM_data: pandas type
        @returns average energy consumption obtained after the start of the experiment
    """
    timelim1,list_with_st_time,energy_consumed = datetime.datetime.strptime(dmt.values[-1][2].split()[3], '%H:%M:%S'),[],[]
    for i in dM_data.values[::-1][1:]:
        if datetime.datetime.strptime(i[2].split()[1], '%H:%M:%S,%f') < timelim1:
            break
        energy_consumed += [float(str(i[6]).replace(',', '.')) * LINE_VOLTAGE * float(str(i[3]).split(":")[2].replace(',', '.'))]
    return sum(energy_consumed)/len(energy_consumed)


def get_data_before_and_after():
    """ Displays the average energy consumption obtained before and after Meter experiment """
    print("\nBefore")
    for dmt,dm in zip(dmt_tup,dm_tup):
        print(check_data_before_meter_start(dmt,dm))
    print("\nAfter")
    for dmt, dm in zip(dmt_tup, dm_tup):
        print(check_data_after_meter_start(dmt, dm))