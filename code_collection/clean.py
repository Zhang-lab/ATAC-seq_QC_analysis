import os
# import time

def init(directory, xpair=None, xsingle=None, xpeak=None, xbed=None):
    if not directory.endswith('/'):
        directory = directory + '/'
    print "Data directory is:", directory
    if not os.path.isdir(directory):
        os.makedirs(directory)
    if xpair is not None:
        pair(directory, xpair)
    if xsingle is not None:
        single(directory, xsingle)
    if xpeak is not None:
        peak(directory, xpeak)
    if xbed is not None:
        bed(directory, xbed)

def pair(directory, input):
    for filename in os.listdir(directory):
        if filename.endswith('pair'):
            os.remove(directory+filename)
    with open(input, 'r') as f:
        for line in f.readlines():
            l = line.split()
            filename = l[0]
            if '_' in filename:
                continue
            with open(directory+filename+'_pair', 'a') as x:
                x.write(line)

def single(directory, input):
    for filename in os.listdir(directory):
        if filename.endswith('single'):
            os.remove(directory+filename)
    with open(input, 'r') as f:
        for line in f.readlines():
            l = line.split()
            filename = l[0]
            if '_' in filename:
                continue
            with open(directory+filename+'_single', 'a') as x:
                x.write(line)

def peak(directory, input):
    for filename in os.listdir(directory):
        if filename.endswith('peak'):
            os.remove(directory+filename)
    with open(input, 'r') as f:
        for line in f:
            l = line.split()
            filename = l[0]
            if '_' in filename:
                continue
            with open(directory+filename+'_peak', 'a') as x:
                x.write(line)

def read_huge(input):
    with open(input, 'r') as f:
        for line in f:
            yield line

def bed(directory, input):
    for filename in os.listdir(directory):
        if filename.endswith('bed'):
            os.remove(directory+filename)
    # with open(input, 'r') as f:
    #     for line in f:
    #         l = line.split()
    #         filename = l[0]
    #         if '_' in filename:
    #             continue
    #         with open(directory+filename+'_bed', 'a') as x:
    #             x.write(line)
    # start = time.time()
    line_of_bed = 0
    with open(input, 'r') as f:
        filename = ''
        filename_content = ''
        #for line_number, line in enumerate(read_huge(input)):
        for line_number, line in enumerate(f):
            l = line.split()[:3]
            if '_' in l[0]:
                continue
            if line_number != 0 and l[0] != filename:
                with open(directory+filename+'_bed', 'w') as x:
                    x.write(filename_content)
                print filename+' writing done...'
                # print 'Time elapsed:', time.time()-start
                filename_content = ''
            filename = l[0].strip()
            filename_content += ' '.join(l)+'\n'
            line_of_bed = line_number
        with open(directory+filename+'_bed', 'w') as x:
            x.write(filename_content)

    with open('config.txt', 'w') as f:
        f.write(str(line_of_bed))
