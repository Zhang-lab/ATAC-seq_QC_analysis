import os
import re
import sys
import clean

OUTPUT = 'peak.txt'
OUTPUT2 = 'reads.txt'

def binary_search(arr, start, end, key, threshold):
    while (start <= end):
        mid = start + (end - start) / 2
        if arr[mid] < key - threshold:
            start = mid + 1
        elif arr[mid] > key + threshold:
            end = mid - 1
        else:
            return mid
    return -1

def clean_file(filename):
    with open(filename, 'w') as f:
        f.write('')

def match(chr, dd_pair, dd_single, output):
    num = 0
    arr_start = [int(x[0]) for x in dd_single]
    arr_end = [int(x[1]) for x in dd_single]
    n = len(dd_single) - 1
    for peak in dd_pair:
        corr_start = binary_search(arr_start, 0, n, int(peak[0]), THRESHOLD)
        corr_end = binary_search(arr_end, 0, n, int(peak[1]), THRESHOLD)
        if corr_start != -1:
            with open(output, 'a') as f:
                f.write(chr+' '+' '.join(peak)+' '+' '.join(dd_single[corr_start])+'\n')
            num += 1
        elif corr_end != -1:
            with open(output, 'a') as f:
                f.write(chr+' '+' '.join(peak)+' '+' '.join(dd_single[corr_end])+'\n')
            num += 1
        else:
            with open('peak1.txt', 'a') as f:
                f.write(chr+' '+' '.join(peak)+'\n')

    arr_start = [int(x[0]) for x in dd_pair]
    arr_end = [int(x[1]) for x in dd_pair]
    n = len(dd_pair) - 1
    for single in dd_single:
        corr_start = binary_search(arr_start, 0, n, int(single[0]), THRESHOLD)
        corr_end = binary_search(arr_end, 0, n, int(single[1]), THRESHOLD)
        if corr_start == -1 and corr_end == -1:
            with open('peak2.txt', 'a') as f:
                f.write(chr+' '+' '.join(single)+'\n')

    return num

def init_match():
    total_len = 0
    total_len2 = 0
    total_match = 0
    clean_file('peak1.txt')
    clean_file('peak2.txt')
    clean_file(OUTPUT)
    for file_pair in m_sort(os.listdir(DIRECTORY)):
        if file_pair.endswith('pair'):
            file_single = file_pair.replace('pair', 'single')
            with open (DIRECTORY+file_pair, 'r') as f:
                dd_pair = []
                for line in f.readlines():
                    l = line.split()[1:4]
                    dd_pair.append(l)
            with open (DIRECTORY+file_single, 'r') as f:
                dd_single = []
                for line in f.readlines():
                    l = line.split()[1:4]
                    dd_single.append(l)
            num = match(re.search(r'(.*?)_',file_pair).group()[:-1],dd_pair, dd_single, OUTPUT)
            total_len += len(dd_pair)
            total_len2 += len(dd_single)
            total_match += num
    # print 1.0 * total_match / total_len
    print total_match, total_len, total_len2

def m_sort(sorted_dir):
    sorted_dir = [x for x in sorted_dir if re.search('chr(.+?)_', x) != None]
    n = len(sorted_dir)
    for i in range(n-2):
        for j in range(i+1, n-1):
            xi = re.search('chr(.+?)_',sorted_dir[i]).group(1)
            xj = re.search('chr(.+?)_',sorted_dir[j]).group(1)
            if xi.isdigit():
                xi = xi.zfill(5)
            if xj.isdigit():
                xj = xj.zfill(5)
            if (xi>xj):
                temp = sorted_dir[i]
                sorted_dir[i] = sorted_dir[j]
                sorted_dir[j] = temp
    return sorted_dir

def count(chr, dd_bed, dd_peak, output):
    index_peak = 0
    n = len(dd_peak)
    reads = [0] * n
    for (x, y) in dd_bed:
        mid = (int(x) + int(y))/2
        if (mid >= int(dd_peak[index_peak][0]) and mid <= int(dd_peak[index_peak][1])):
            reads[index_peak] += 1
        elif (index_peak<n-1) and (mid >= int(dd_peak[index_peak+1][0]) and mid <= int(dd_peak[index_peak+1][1])):
            index_peak += 1
            reads[index_peak] += 1
        elif (index_peak<n-2) and (mid >= int(dd_peak[index_peak+2][0]) and mid <= int(dd_peak[index_peak+2][1])):
            index_peak += 2
            reads[index_peak] += 1
    with open(output, 'a') as f:
        for index_peak, x in enumerate(reads):
            if x != 0 and index_peak < n:
                y = 1.0*x*10**9/LINE_OF_BED/(int(dd_peak[index_peak][1])-int(dd_peak[index_peak][0]))
                y = '%.4f' %y
                f.write(chr+' '+' '.join(dd_peak[index_peak])+' '+str(x)+' '+y+' '+'\n')
            else:
                if index_peak < n:
                    f.write(chr+' '+' '.join(dd_peak[index_peak])+' '+str(x)+' '+'0'+' '+'\n')
    
def init_count():
    clean_file(OUTPUT2)
    for file_bed in m_sort(os.listdir(DIRECTORY)):
        if file_bed.endswith('bed'):
            file_peak = file_bed.replace('bed', 'peak')
            with open (DIRECTORY+file_bed, 'r') as f:
                dd_bed = []
                for line in f.readlines():
                    l = line.split()[1:3]
                    dd_bed.append(l)
            if not os.path.isfile(DIRECTORY+file_peak):
                continue
            with open (DIRECTORY+file_peak, 'r') as f:
                dd_peak = []
                for line in f.readlines():
                    try:
                        l = line.split()[1:6]
                        l[0] = str(min(int(l[0]), int(l[3])))
                        l[1] = str(max(int(l[1]), int(l[4])))
                        ll = l[:3]
                    except:
                        ll = line.split()[1:4]
                    dd_peak.append(ll)
            print file_bed+" reading done..."
            count(re.search(r'(.*?)_',file_bed).group()[:-1], dd_bed, dd_peak, OUTPUT2)
            #except:

if len(sys.argv) == 5:
    THRESHOLD = int(sys.argv[1])
    DIRECTORY = sys.argv[2]
    clean.init(directory=DIRECTORY, xpair=sys.argv[3], xsingle=sys.argv[4])
    init_match()

if len(sys.argv) == 4:
    DIRECTORY = sys.argv[1]
    clean.init(directory=DIRECTORY, xpeak=sys.argv[3], xbed=sys.argv[2])
    with open('config.txt', 'r') as f:
        LINE_OF_BED = int(f.read().strip())
    init_count()
