import pandas as pd
import datetime

############# METERS
#### DAY 1

# Read the CSV file
dm = pd.read_csv('dataMeters.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt = pd.read_csv('dataMetersTime.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])

dm2 = pd.read_csv('dataMeters2.csv', sep=';', usecols=['Reading', 'Sample', 'Start Time', 'Duration', 'Max Time', 'Max', 'Average', 'Min Time', 'Min', 'Description', 'Stop Time'])
dmt2 = pd.read_csv('dataMetersTime2.csv', sep=',', usecols=['Configurations', 'Start Time', 'Stop Time'])


#print(dm.loc[range(5),:])
#print(dmt.loc[range(5),:])

# starts_dmt = [i.split(" ")[4] for i in dmt["Start Time"]]
# ends_dmt = starts_dmt[1:]+[(dmt["Stop Time"][-1:].values[0]).split(" ")[4]]
# starts_dm = [i.split(" ")[1] for i in dm["Start Time"]]
# ends_dm = starts_dm[1:]+[(dm["Stop Time"][-1:].values[0]).split(" ")[1]]
# startsends_dmt= zip(starts_dmt,ends_dmt)
# startsends_dm= zip(starts_dm,ends_dm)

## optimization

def get_zip_startend_time(source,ind):
    starts = [i.split(" ")[ind] for i in source["Start Time"]]
    ends = starts[1:] + [(source["Stop Time"][-1:].values[0]).split(" ")[ind]]
    return zip(starts,ends)

startsends_dmt= get_zip_startend_time(dmt,4)
startsends_dm = get_zip_startend_time(dm,1)

# print(list(startsends_dmt))
# print(list(startsends_dm))

## Work with formulas

# power_in_watts = [float(i.replace(",","."))*237.2 for i in dm['Average'] if str(i) !="nan"]
# duration_in_seconds = [(datetime.datetime.strptime(i.replace(",","."), '%H:%M:%S.%f').second +datetime.datetime.strptime(i.replace(",","."), '%H:%M:%S.%f').microsecond)/1e6  for i in dm['Duration'] if str(i) !="nan"]
# energy_in_joules = [p*d for (p,d) in zip(power_in_watts,duration_in_seconds)]
#
#print(sum(energy_in_joules))


#('08:50:14', '08:50:20'), ('08:50:20', '08:50:24'), ('08:50:24', '08:50:29'), ('08:50:29', '08:50:33'), ('08:50:33', '08:50:37'), ('08:50:37', '08:50:41'), ('08:50:41', '08:50:46'),
#('8:49:12,1', '8:49:12,9'), ('8:49:12,9', '8:49:14,1'), ('8:49:14,1', '8:49:14,7'), ('8:49:14,7', '8:49:15,3'), ('8:49:15,3', '8:49:16,1'), ('8:49:16,1', '8:49:20,4'), ('8:49:20,4', '8:49:25,4')

## attempt 0

# def is_time_within_frame(time_b, frame_a):
#     start_time, end_time = frame_a
#     return start_time <= time_b <= end_time
# tr=0
# for time_b in startsends_dm:
#     for frame_a in startsends_dmt:
#         if is_time_within_frame(time_b[0], frame_a) or is_time_within_frame(time_b[1], frame_a):
#             print(tr, f"{time_b} falls within {frame_a}")
#             tr+=1
#('9:05:42,5', '9:05:44,6')
# a_dt = [(datetime.datetime.strptime(start, '%H:%M:%S'), datetime.datetime.strptime(end, '%H:%M:%S')) for start, end in startsends_dmt]
#
# # Convert time strings in b to datetime objects
# b_dt = [(datetime.datetime.strptime(start, '%H:%M:%S,%f'),datetime.datetime.strptime(end, '%H:%M:%S,%f')) for start, end in startsends_dm]
# count =0

# # Iterate over each pair of datetime objects in a
# for ca, (start_a, end_a) in enumerate(a_dt):
#     # Check if any datetime object in b falls within the time frame
#     for cb,(start_b, end_b) in enumerate(b_dt):
#         if start_a <= start_b and  end_b<= end_a:
#             # If datetime object from b falls within the time frame, do something with it
#             print(count," -- ",ca, " __ ",str(start_b).split(" ")[1]," >< ",str(end_b).split(" ")[1]," -- ",cb, " __ ",str(start_a).split(" ")[1]," >< ",str(end_a).split(" ")[1])
#             count+=1
#             break

#### DAY 2

## attempt 1

# def find_matching_rows(a, b):
#     matches = {}
#     for start_a, end_a in a:
#         best_match = None
#         for start_b, end_b in b:
#             if start_a <= start_b <= end_a:
#                 if not best_match or start_b < best_match[0]:
#                     best_match = (start_b, end_b)
#         if best_match:
#             matches[(start_a, end_a)] = best_match
#     return matches

# print(len(b_dt) , len((matches.values())))
# print(list((matches.keys()))[0])
# print(list((matches.values()))[0])
# print(" xxxxx ")
# print(list((matches.keys()))[1])
# print(list((matches.values()))[1])
# print(" xxxxx ")
# print(list((matches.keys()))[2])
# print(list((matches.values()))[2])
# print(type(list(matches.keys())[0][0]))
# print(list(matches.keys())[0][0])

## attempt 2

# def get_matching_frames(a, b):
#     b_dict = {t: row for row in b for t in row}
#     matching_frames = {}
#     ids = []
#     for start, end in a:
#         if (start, end) in b_dict:
#             matching_frames[(start, end)] = b_dict[(start, end)]
#         else:
#             closest_frame_start = min(t for t in b_dict if t >= start)
#             matching_frames[(start, end)] = b_dict[closest_frame_start]
#             ids.append(list(b_dict.keys()).index(closest_frame_start))
#     return (matching_frames,ids)

## attempt 3

# def get_matching_frames2(a, b):
#     matching_frames,ids,c = [],[],0
#     for (start_a,end_a) in a:
#         lenOfIds = len(ids)
#         for (start_b, end_b) in b:
#             c+=1
#             if start_a < start_b and end_b < end_a:
#                 matching_frames.append((start_b,end_b))
#                 b = b[b.index((start_b, end_b)):]
#                 ids.append(c)
#                 break
#         if lenOfIds == len(ids):
#             closest_frame_start = min(t for t in b[:7] if t[0] >= start_a)
#             matching_frames.append(closest_frame_start)
#             ids.append(c)
#     return (matching_frames,ids)

## attempt 4

def get_matching_frames3(a, b):
    matching_frames,ids = [],[]
    mega_b = b
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

# date_format = "%Y-%m-%d %H:%M:%S"
#
# date_string = list(matches.keys())[0][0].strftime(date_format)
# print(date_string)
# a = [('09:50:43', '09:50:46'), ('09:50:46', '09:50:49'), ('09:50:49', '09:50:53'), ('09:50:53', '09:50:56'), ('09:50:56', '09:50:59'), ('09:50:59', '09:51:03')]
# b = [('9:50:40,1', '9:50:40,4'), ('9:50:40,4', '9:50:45,2'), ('9:50:45,2', '9:50:45,3'), ('9:50:45,3', '9:50:45,9'), ('9:50:45,9', '9:50:47,4'), ('9:50:47,4', '9:50:49,3'), ('9:50:49,3', '9:50:50,4'), ('9:50:50,4', '9:50:52,6'), ('9:50:52,6', '9:50:54,1'), ('9:50:54,1', '9:50:55,3'), ('9:50:55,3', '9:50:55,4'), ('9:50:55,4', '9:51:00,1'), ('9:51:00,1', '9:51:00,3'), ('9:51:00,3', '9:51:02,6')]

a_dt = [(datetime.datetime.strptime(start, '%H:%M:%S'), datetime.datetime.strptime(end, '%H:%M:%S')) for start, end in startsends_dmt]
b_dt = [(datetime.datetime.strptime(start, '%H:%M:%S,%f'),datetime.datetime.strptime(end, '%H:%M:%S,%f')) for start, end in startsends_dm]

(matches,ids) = get_matching_frames3(a_dt, b_dt)

## testing

#print(dm.loc[ids[0]])
# print(len(matches))
# print(len(set(matches)))
#
# print(ids)
# print(len(ids))
# print(len(set(ids)))
#
# povtor = set([i for i in ids if ids.count(i) >1 ])
# print(povtor)
# print(len(povtor))
# # print(ids)
# # for c,i in enumerate(ids):
# #     print(dm.loc[i,:])
#

## results

power_in_watts = [float(i.replace(",","."))*237.2 for c,i in enumerate(dm['Average']) if str(i) !="nan" and c in ids]
duration_in_seconds = [(datetime.datetime.strptime(i, '%H:%M:%S,%f').second +datetime.datetime.strptime(i, '%H:%M:%S,%f').microsecond)/1e6  for c,i in enumerate(dm['Duration']) if str(i) !="nan" and c in ids]

energy_in_joules = [p*d*10 for (p,d) in zip(power_in_watts,duration_in_seconds)]
count=1

for count in range(int(len(energy_in_joules)/33)):
    print(count+1," -- ", sum(energy_in_joules[:33])/33)
    energy_in_joules= energy_in_joules[33:]


## data view

#k = 0
# for c,(i,j) in enumerate(zip(a_dt,matches)):
#     if i[1]<j[0]:
#         k+=1
#         print(c+1, "Start: ",str(i[0]).split(" ")[1], " <= ", str(j[0]).split(" ")[1], " -- ", str(j[1]).split(" ")[1] ," <= ", str(i[1]).split(" ")[1], ": Stop")
# print(k)


#
# print(len(mass_with_ress_meters))
# #for c,i in enumerate(list(mass_fir.values())+list(mass_sec.values())): mass_with_ress_meters[c+1] = i
#
#
# print(len(mass_with_ress_rapl))
# # print(mass_with_ress_meters)
# # print(mass_with_ress_rapl)
