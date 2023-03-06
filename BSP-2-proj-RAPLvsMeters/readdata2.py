import pandas as pd
import matplotlib.pyplot as plt
import datetime

mass_threads= [1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64]

dm = pd.read_csv('dataMeters.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt = pd.read_csv('dataMetersTime.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])

dm2 = pd.read_csv('dataMeters2.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt2 = pd.read_csv('dataMetersTime2.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])

dr = pd.read_csv('dataRAPL.csv', sep=',', usecols=['version', 'n-workers', 'n-secondary-workers', 'n-reservations', 'n-relations', 'n-queries', 'password-work-factor', 'Energy Cores', 'Energy Ram', 'Energy Gpu', 'Energy Pkg', 'Time (perf) (s)', 'Time (exec) (s)', 'Ratio time (%)', 'Power_Cores', 'Power_Ram', 'Power_Gpu', 'Power_Pkg'])

def get_zip_startend_time(source,ind):
    starts = [i.split(" ")[ind] for i in source["Start Time"]]
    ends = starts[1:] + [(source["Stop Time"][-1:].values[0]).split(" ")[ind]]
    return zip(starts,ends)

def get_data_for_15(a, b,source_dm):
    id_dict, energy_ress,mega_b = 1, {i: [] for i in mass_threads},b
    for (start_a, end_a) in a:
        #print(">>",id_dict, " <<")
        #print("    ", start_a," --- ", end_a)
        for (start_b, end_b) in b:
            if start_a <= start_b and end_b <= end_a:
                energy_ress[id_dict] +=[float(str(source_dm["Average"][mega_b.index((start_b, end_b))]).replace(',','.'))*237.2 * float(str(source_dm["Duration"][mega_b.index((start_b, end_b))]).split(":")[2].replace(',','.'))]
                #print(float(str(source_dm["Duration"][b.index((start_b, end_b))]).split(":")[2].replace(',','.')))
                #print("    "*2,start_b, " -- ", end_b, " ~ ", float(str(source_dm["Average"][b.index((start_b, end_b))]).replace(',','.')), " * ",237.2, " * ",  float(str(source_dm["Duration"][b.index((start_b, end_b))]).split(":")[2].replace(',','.')))
            if start_b>start_a and end_b > end_a:
                b = b[b.index((start_b, end_b))-1:]
                break
        id_dict+=1 if id_dict ==1 else 2
        if id_dict == 66: id_dict = 1
    return [sum(energy_ress[i])/len(energy_ress[i]) for i in mass_threads]

def get_time_in_time(dm,dmt):
    a_dt_1 = [(datetime.datetime.strptime(start, '%H:%M:%S'), datetime.datetime.strptime(end, '%H:%M:%S')) for start, end in get_zip_startend_time(dmt,4)]
    b_dt_1 = [(datetime.datetime.strptime(start, '%H:%M:%S,%f'),datetime.datetime.strptime(end, '%H:%M:%S,%f')) for start, end in get_zip_startend_time(dm,1)]
    return (a_dt_1,b_dt_1)


def get_res_data_rapl(dr):
    return [i + j + k + t for (i, j, k, t) in zip(dr['Energy Cores'], dr['Energy Ram'], dr['Energy Gpu'], dr['Energy Pkg'])]

def calc_avar(mass,source):
    for j in range(33):mass.append(sum(list(source[j::33])) /(len(source)/ 33))

(a_dt_1,b_dt_1) = get_time_in_time(dm,dmt)
(a_dt_2,b_dt_2) = get_time_in_time(dm2,dmt2)

mega_mass_with_mega_data ={(i+1):(get_data_for_15(a_dt_1[495*i:495+495*i], b_dt_1, dm) ,get_data_for_15(a_dt_2[495 * i:495 + 495 * i], b_dt_2, dm2))  for i in range(5)}
newMMWMD = {}


for j in range(1,6):
    sumpi = [(mega_mass_with_mega_data[j][0][i] + mega_mass_with_mega_data[j][1][i])/2 for i in range(33)]
    newMMWMD[j] =sumpi

mass_with_ress = get_res_data_rapl(dr)
mass_data_rapl_or =[]
mass_data_rapl_tx1 =[]
mass_data_rapl_tx2 =[]
mass_data_rapl_tx3 =[]
mass_data_rapl_tx4 =[]

calc_avar(mass_data_rapl_or,mass_with_ress[0:990])
calc_avar(mass_data_rapl_tx1,mass_with_ress[990:1980])
calc_avar(mass_data_rapl_tx2,mass_with_ress[1980:2970])
calc_avar(mass_data_rapl_tx3,mass_with_ress[2970:3960])
calc_avar(mass_data_rapl_tx4,mass_with_ress[3960:4950])

yl = [1]+list(range(2,66,2))

figure, axis = plt.subplots(2)
#

# axis[0].plot(yl, mega_mass_with_mega_data[1][0],"b",label= "origin")
# axis[0].plot(yl, mega_mass_with_mega_data[2][0],"r",label= "txact")
# axis[0].plot(yl, mega_mass_with_mega_data[3][0],"r",label= "txact")
# axis[0].plot(yl, mega_mass_with_mega_data[4][0],"r",label= "txact")
# axis[0].plot(yl, mega_mass_with_mega_data[5][0],"r",label= "txact")
# axis[0].set_title("Meters data (origin blue) energy consumtion 15")
# axis[0].set_xlabel('1-64 threads')
# axis[0].set_ylabel('energy consumption (J)')
#

# axis[1].plot(yl, mega_mass_with_mega_data[1][1],"b",label= "origin")
# axis[1].plot(yl, mega_mass_with_mega_data[2][1],"r",label= "txact")
# axis[1].plot(yl, mega_mass_with_mega_data[3][1],"r",label= "txact")
# axis[1].plot(yl, mega_mass_with_mega_data[4][1],"r",label= "txact")
# axis[1].plot(yl, mega_mass_with_mega_data[5][1],"r",label= "txact")
# axis[1].set_title("Meters data (origin blue) energy consumtion 15")
# axis[1].set_xlabel('1-64 threads')
# axis[1].set_ylabel('energy consumption (J)')

# For Cosine Function

axis[0].plot(yl, newMMWMD[1],"b",label= "origin")
axis[0].plot(yl, newMMWMD[2],"r",label= "txact")
axis[0].plot(yl, newMMWMD[3],"r",label= "txact")
axis[0].plot(yl, newMMWMD[4],"r",label= "txact")
axis[0].plot(yl, newMMWMD[5],"r",label= "txact")
axis[0].set_title("Meters data (origin blue) energy consumtion 30")
#axis[0].set_xlabel('1-64 threads')
axis[0].set_ylabel('energy consumption (J)')
axis[0].legend()

axis[1].plot(yl, mass_data_rapl_or,"m",label= "origin")
axis[1].plot(yl, mass_data_rapl_tx1,"c",label= "txact")
axis[1].plot(yl, mass_data_rapl_tx2,"c",label= "txact")
axis[1].plot(yl, mass_data_rapl_tx3,"c",label= "txact")
axis[1].plot(yl, mass_data_rapl_tx4,"c",label= "txact")
#axis[1].set_title("RAPL data (origin magenta) energy consumtion 30")
axis[1].set_xlabel('RAPL data (origin magenta) energy consumtion 30')
axis[1].set_ylabel('energy consumption (J)')
# Combine all the operations and display
axis[1].legend()
plt.show()




