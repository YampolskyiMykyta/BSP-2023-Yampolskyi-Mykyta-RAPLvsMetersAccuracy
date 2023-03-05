import pandas as pd
import matplotlib.pyplot as plt
import datetime

def get_zip_startend_time(source,ind):
    starts = [i.split(" ")[ind] for i in source["Start Time"]]
    ends = starts[1:] + [(source["Stop Time"][-1:].values[0]).split(" ")[ind]]
    return zip(starts,ends)

def get_matching_frames3(a, b):
    matching_frames,ids,mega_b = [],[],b
    for (start_a, end_a) in a:
        bip=0
        for (start_b, end_b) in b:
            bip+=1
            if start_a < start_b and end_b < end_a:
                matching_frames.append((start_b,end_b))
                b = b[b.index((start_b, end_b)):]
                ids.append(mega_b.index((start_b, end_b)))
                break
            if bip == len(b):
                for i in b[0:10]:
                    if mega_b.index(i) not in ids:
                        matching_frames.append(i)
                        ids.append(mega_b.index(i))
                        break
    return (matching_frames,ids)

def get_res_data_meters(dm, dmt):
    a_dt = [(datetime.datetime.strptime(start, '%H:%M:%S'), datetime.datetime.strptime(end, '%H:%M:%S')) for start, end in get_zip_startend_time(dmt,4)]
    b_dt = [(datetime.datetime.strptime(start, '%H:%M:%S,%f'),datetime.datetime.strptime(end, '%H:%M:%S,%f')) for start, end in get_zip_startend_time(dm,1)]

    (matches,ids) = get_matching_frames3(a_dt, b_dt)

    #matches == good
    # for c,i in enumerate(matches):
    #     print(c+1 , "  ",str(a_dt[c][0]).split(" ")[1]," -- ", " < ",str(i[0]).split(" ")[1],"  -  ", str(i[1]).split(" ")[1], " < ",str(a_dt[c][1]).split(" ")[1])

    # ids == good
    # print(len(ids))
    # for c,i in enumerate(matches):
    #     print(c+1 , "  ", str(dm["Start Time"][ids[c]]).split(" ")[1] , " = ",str(i[0]).split(" ")[1],"  -  ", str(i[1]).split(" ")[1], " = ", str(dm["Stop Time"][ids[c]]).split(" ")[1])

    power_in_watts = [float(f"0.{i.replace(',','')}")*237.2 for c,i in enumerate(dm['Average']) if str(i) !="nan" and c in ids]
    #power_in_watts = [float(i.replace(",","."))*237.2 for c,i in enumerate(dm['Average']) if str(i) !="nan" and c in ids]
    duration_in_seconds = [(datetime.datetime.strptime(i, '%H:%M:%S,%f').second +datetime.datetime.strptime(i, '%H:%M:%S,%f').microsecond)/1e6  for c,i in enumerate(dm['Duration']) if str(i) !="nan" and c in ids]

    return [p*d for (p,d) in zip(power_in_watts,duration_in_seconds)]

def get_res_data_rapl(dr):
    return [i + j + k + t for (i, j, k, t) in zip(dr['Energy Cores'], dr['Energy Ram'], dr['Energy Gpu'], dr['Energy Pkg'])]

def calc_avar(mass,source):
    for j in range(33):mass.append(sum(list(source[j::33])) /(len(source)/ 33))

def show_and_get_differences(rapl,meters,toshow):
    diffs =[abs(r-m) for (r,m) in zip(rapl,meters) ]
    if(toshow):
        print()
        print(" NT |", " " * 7, "RAPL", " " * 8, "|", " " * 8, "METERS"," "*6,"|"," "*4,"DIFFERENCES")
        print("-" * 76)
        for i in range(len(meters)):
            if ((i + 1) <= 2):
                print('{:>2}'.format(i + 1), " | ", '{:>20}'.format(rapl[i]), "|", '{:>20}'.format(meters[i]), f"  |   {diffs[i]}" )
            else:
                print('{:>2}'.format(((i + 1) * 2) - 2), " | ", '{:>20}'.format(rapl[i]), "|", '{:>20}'.format(meters[i]), f"  |   {diffs[i]}")
        print()
        return diffs

dm = pd.read_csv('dataMeters.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt = pd.read_csv('dataMetersTime.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])

dm2 = pd.read_csv('dataMeters2.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt2 = pd.read_csv('dataMetersTime2.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])

dr = pd.read_csv('dataRAPL.csv', sep=',', usecols=['version', 'n-workers', 'n-secondary-workers', 'n-reservations', 'n-relations', 'n-queries', 'password-work-factor', 'Energy Cores', 'Energy Ram', 'Energy Gpu', 'Energy Pkg', 'Time (perf) (s)', 'Time (exec) (s)', 'Ratio time (%)', 'Power_Cores', 'Power_Ram', 'Power_Gpu', 'Power_Pkg'])


mass_with_ress_meters = get_res_data_meters(dm, dmt[:495]) + get_res_data_meters(dm2, dmt2[:495])
mass_with_ress_rapl = get_res_data_rapl(dr[:990])
#
# print(len(mass_with_ress_meters))
# print(mass_with_ress_rapl[::33])
#
# # print(min(mass_with_ress_rapl[::33]))
# # print(max(mass_with_ress_rapl[::33]))
# # print((mass_with_ress_rapl[::33]))
#
# mass_avar_meters =[]
# calc_avar(mass_avar_meters, mass_with_ress_meters)
#
# mass_avar_rapl =[]
# calc_avar(mass_avar_rapl, mass_with_ress_rapl)
#
# diff = show_and_get_differences(mass_avar_rapl,mass_avar_meters,True)
# print(max(diff))
#
# yl = [1]+list(range(2,66,2))

# plt.plot(yl, mass_avar_rapl,'ro-')
# plt.plot(yl, mass_avar_meters,'bo-')
# plt.title("RAPL (red) vs Meters (blue) energy consumtion")
# plt.xlabel('1-64 threads')
# plt.ylabel('energy consumption (J)')
# plt.show()


