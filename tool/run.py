# run.py

#/* --------------------------------------- */
#/* This script is written by Sunghye Park  */
#/* 2022.03.16                              */
#/*                                         */
#/* shpark96@postech.ac.kr                  */
#/* --------------------------------------- */

import os
import signal
import sys
import shlex
import subprocess as sp
import argparse as ap
import multiprocessing as mp
from datetime import datetime
from time import time as timer, sleep

# The Qmap mapper
RESULT_QMAP = (
    ("alu_bdd_288"         , 109     , 15      , 13  , 254    , 1.77362 ),
    ("alu_v0_27"           , 34      , 3       , 4   , 106    , 4.13529 ),
    ("benstein_vazirani"   , 1       , 0       , 0   , 36     , 0.01135 ),
    ("4gt12_v1_89"         , 259     , 51      , 3   , 690    , 26.0342 ),
    ("4gt4_v0_72"          , 273     , 52      , 2   , 788    , 4.628   ),
    ("4mod5_bdd_287"       , 69      , 10      , 4   , 226    , 19.5225 ),
    ("cm42a_207"           , 2301    , 494     , 24  , 5711   , 535.889 ),
    ("cnt3_5_180"          , 641     , 142     , 0   , 1235   , 25.6028 ),
    ("cuccaroAdder_1b"     , 23      , 0       , 3   , 90     , 0.28906 ),
    ("cuccaroMultiply"     , 64      , 6       , 7   , 216    , 2.09601 ),
    ("decod24_bdd_294"     , 83      , 15      , 3   , 201    , 1.46449 ),
    ("decod24_enable_126"  , 434     , 95      , 0   , 1150   , 28.7617 ),
    ("graycode6_47"        , 5       , 0       , 0   , 16     , 5.8601  ),
    ("ham3_102"            , 17      , 2       , 0   , 60     , 0.22297 ),
    ("miller_11"           , 39      , 0       , 8   , 139    , 0.27471 ),
    ("mini_alu_167"        , 294     , 56      , 0   , 818    , 28.5271 ),
    ("mod5adder_127"       , 598     , 109     , 16  , 1617   , 74,544  ),
    ("mod8_10_177"         , 518     , 106     , 2   , 1436   , 2.53567 ),
    ("one_two_three"       , 70      , 10      , 4   , 215    , 6.41456 ),
    ("rd32_v0_66"          , 31      , 1       , 6   , 104    , 1.71454 ),
    ("rd53_311"            , 375     , 81      , 4   , 856    , 0.67113 ),
    ("rd73_140"            , 292     , 52      , 16  , 662    , 20.3441 ),
    ("rd84_142"            , 448     , 98      , 0   , 861    , 21.1735 ),
    ("sf_274"              , 818     , 104     , 85  , 2150   , 41.2879 ),
    ("shor_15"             , 4924    , 982     , 95  , 11218  , 14.6928 ),
    ("sqrt8_260"           , 4216    , 956     , 17  , 10016  , 13.2009 ),
    ("squar5_261"          , 2630    , 585     , 3   , 6468   , 7.76352 ),
    ("sym6_145"            , 4787    , 970     , 88  , 12869  , 16.7757 ),
    ("sym9_146"            , 431     , 91      , 5   , 980    , 21.6935 ),
    ("sys6_v0_111"         , 267     , 49      , 11  , 4389   , 21.1608 ),
    ("vbeAdder_2b"         , 80      , 6       , 10  , 215    , 0.1938  ),
    ("wim_266"             , 1254    , 265     , 16  , 3190   , 13.5583 ),
    ("xor5_254"            , 8       , 1       , 0   , 18     , 29.2384 ),
    ("z4_268"              , 4088    , 887     , 42  , 9703   , 869.445 ),
    ("adr4_197"            , 4685    , 1021    , 62  , 11070  , 10.924  ),
    ("9symml_195"          , 49154   , 11282   , 38  , 116097 , 2332.7  ),
    ("clip_206"            , 49090   , 11268   , 257 , 111277 , 2587.99 ),
    ("cm152a_212"          , 1591    , 353     , 0   , 4086   , 3.53968 ),
    ("cm85a_209"           , 16224   , 3716    , 45  , 37827  , 389.036 ),
    ("co14_215"            , 27787   , 6615    , 51  , 57972  , 1044.04 ),
    ("cycle10_2_110"       , 8471    , 1899    , 63  , 20449  , 106.071 ),
    ("dc1_220"             , 2444    , 481     , 84  , 5978   , 8.51783 ),
    ("dc2_222"             , 13520   , 3077    , 79  , 31796  , 268.637 ),
    ("dist_223"            , 54717   , 12599   , 148 , 124030 , 3550.7  ),
    ("ham15_107"           , 12030   , 2704    , 30  , 28896  , 193.48  ),
    ("life_238"            , 31920   , 7324    , 74  , 75422  , 1364.49 ),
    ("max46_240"           , 37565   , 8217    , 535 , 89167  , 1840.89 ),
    ("mini_alu_305"        , 242     , 41      , 21  , 518    , 0.40925 ),
    ("misex1_241"          , 6588    , 1480    , 24  , 15884  , 53.4152 ),
    ("pm1_249"             , 2331    , 504     , 24  , 5627   , 6.86838 ),
    ("radd_250"            , 4363    , 948     , 57  , 10798  , 24.0061 ),
    ("root_255"            , 24520   , 5627    , 73  , 55952  , 844.114 ),
    ("sqn_258"             , 13815   , 2984    , 202 , 33000  , 270.981 ),
    ("square_root_7"       , 9845    , 2184    , 102 , 23201  , 46862.9 ),
    ("sym10_262"           , 92326   , 20986   , 642 , 215222 , 9083.41 ),
    ("sym9_148"            , 28462   , 6182    , 254 , 70748  , 800.226 ),
                )
# define paths
dirpos = "../example"
binaryName = "./MCQA"
outpos = "../output"
logpos = "../log"

#
curTime = datetime.now().strftime('%m_%d_%H_%M')
benchlist_dir = []
benchlist = []
benchDir = sorted(os.listdir(dirpos))

# result include tuple
benchresult = []

#### function define ####
def ExecuteCommand( curBench, result, alpha, beta, time_limit):
    inputPath = "%s/%s" % (dirpos, curBench)
    outputPath = "%s/%s" % (outpos, curBench)
    logPath = "%s/%s.%s.log" % (logpos, curBench, curTime)
    exeStr = "%s %s %s %s %s | tee %s" % (binaryName, inputPath, outputPath, alpha, beta, logPath)
    print( exeStr )
    run = RunCommand(exeStr, time_limit)
    if(run!=None):
        ourResult  = ParseLog(logPath)
        TCADResult = FindTCAD(curBench[:-5])
        result.append( [(curBench[:-5], alpha, beta), TCADResult, ourResult] )
        print('Successfully run!')
    else:
        print('Time out!')

def RunCommand(exeStr, time_limit):
    try:
        proc = sp.Popen(exeStr, shell=True)
        return proc.communicate(timeout = time_limit)
    except sp.TimeoutExpired:
        print('Time out!')
        pid = proc.pid
        proc.terminate()
        proc.kill()
        return None

def FindTCAD( name ):
    cz, swap, mov, latency, time = 0, 0, 0, 0, 0
    for rslt in RESULT_QMAP:
        curBench = rslt[0]
        if curBench == name.replace('-', '_'):
            cz      = rslt[1]
            swap    = rslt[2]
            mov     = rslt[3]
            latency = rslt[4]
            time    = rslt[5]
    return (cz, swap, mov, latency, time)

def ParseLog( logPath ):
    cnot, cz, swap, mov, bridge, latency, time = 0, 0, 0, 0, 0, 0, 0
    f = open(logPath)
    for line in f:
        if line[0:16] == "# add_cnot_num: " :
            cnot  = int(line[16:])
        if line[0:14] == "# add_cz_num: " :
            cz  = int(line[14:])
        if line[0:16] == "# add_swap_num: " :
            swap  = int(line[16:])
        if line[0:15] == "# add_mov_num: " :
            mov   = int(line[15:])
        if line[0:18] == "# add_bridge_num: " :
            bridge   = int(line[17:])
        if line[0:11] == "# latency: " :
            latency  = int(line[10:])
        if line[0:11] == "END Program" :
            time = line.rstrip('\n').split(' ')[-1]
    f.close()
    cnot_num = [ cnot, cz, swap, mov, bridge, latency, time ]
    return cnot_num

def getdirList( ):
    for curBench in benchDir:
        if curBench[-5:] == ".qasm":
            benchlist_dir.append(curBench)

##########################
getdirList()
if len(sys.argv) <= 1 or sys.argv[1] == "help":
    print("usage: python run.py <arg>")
    print("")
    print("list of <arg>")
    print("     help        : print how to use run.py")
    print("     list        : print file number of bench directory")
    print("     all         : execute all bench file in directory")
    print("     #           : execute one bench file(or more) in directory")
    print("                   to look file number do list")
    print("     <benchname> : execute one bench file of benchname")
    sys.exit(1)
elif sys.argv[1] == "list":
    for index, curBench in enumerate(benchlist_dir):
        print( str(index) + " : " + curBench )
    sys.exit(1)
elif sys.argv[1].isdigit() and len(sys.argv) == 2:
    benchNum = int(sys.argv[1])
    try:
        benchlist_dir[benchNum]
    except:
        print("ERROR: NO EXISTING BENCH FILE")
        print("Please do python run.py list")
        sys.exit(1)
    ExecuteCommand( benchlist_dir[benchNum], benchresult, 1.0, 3.0, 50) #alpha, beta, time_limit
elif sys.argv[1] == "all" or len(sys.argv) >= 3:
        if len(sys.argv) >= 3:
            for argnum in range(1, len(sys.argv)):
                num = int(sys.argv[argnum])
                try:
                    benchlist_dir[num]
                except:
                    print("ERROR: NO EXISTING BENCH FILE")
                    print("Please do python run.py list")
                    sys.exit(1)
                benchlist.append(benchlist_dir[num])
        elif sys.argv[1] == "all":
            benchlist = benchlist_dir

        if __name__ == '__main__':
            procs = []
            manager = mp.Manager()
            results = manager.list()
            for index, curBench in enumerate(benchlist):
                proc = mp.Process(target=ExecuteCommand,
                        args=(curBench, results, 0.9, 3.0, 50)) #Best parameter: alpha:0.9, beta:3.0
                procs.append(proc)
                proc.start()
            for proc in procs:
                proc.join()
        benchresult = results[:]
else:
    print("ERROR: NO ACCEPTABLE ARGUMENT")
    sys.exit(1)

benchresult = sorted(benchresult, key=lambda row: row[1][0])
avg_cz  = 0
avg_lty = 0.0
avg_rt  = 0.0
print(" ")
print(" ")
#print("----------------------------------------- Result ------------------------------------------")
#print("                    |                                Ours                                 |")
#print("   bench_name       |  #CNOT  |  #CZ    |  #SWAP  |  #MOV   | #BRIDGE | latency | RunTime |")
#print("-------------------------------------------------------------------------------------------")
#for benchResult in benchresult:
#    s1 = '%-19s' % benchResult[0][0]
#    s2 = '%7s' % benchResult[1][0]
#    s3 = '%7s' % benchResult[1][1]
#    s4 = '%7s' % benchResult[1][2]
#    s5 = '%7s' % benchResult[1][3]
#    s6 = '%7s' % benchResult[1][4]
#    s7 = '%7s' % benchResult[1][5]
#    s8 = '%7s' % benchResult[1][6]
#    try:
#        print( s1, "|", s2, "|", s3, "|", s4, "|", s5, "|", s6, "|", s7, "|", s8, "|")
#    except ZeroDivisionError as e:
#        print( s1, e)
#
#print("-------------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------- Result -------------------------------------------------------------------------------")
print("                    ||                     TCAD                        ||                           Ours                            ||          Comparison         |")
print("   bench_name       ||   #CZ   |  #SWAP  |  #MOV   | latency | RunTime ||  #CZ    |  #SWAP  |  #MOV   | #BRIDGE | latency | RunTime ||  #CZ    | latency | RunTime |")
print("--------------------------------------------------------------------------------------------------------------------------------------------------------------------")
for benchResult in benchresult:
    s1 = '%-19s' % benchResult[0][0]
    t1 = '%7s' % benchResult[1][0]
    t2 = '%7s' % benchResult[1][1]
    t3 = '%7s' % benchResult[1][2]
    t4 = '%7s' % benchResult[1][3]
    t5 = '%7s' % benchResult[1][4]
    #
    s3 = '%7s' % benchResult[2][1]
    s4 = '%7s' % benchResult[2][2]
    s5 = '%7s' % benchResult[2][3]
    s6 = '%7s' % benchResult[2][4]
    s7 = '%7s' % benchResult[2][5]
    s8 = '%7s' % benchResult[2][6]
    #
    c_cz  = int(s3)/int(t1)
    c_lty = float(s7)/float(t4)
    c_rt  = float(s8)/float(t5)
    avg_cz  += c_cz
    avg_lty += c_lty
    avg_rt  += c_rt

    try:
        print( s1, "||", 
                t1, "|", t2, "|", t3, "|", t4, "|", t5, "||",
                s3, "|", s4, "|", s5, "|", s6, "|", s7, "|", s8, "||",
                "%7.2f | %7.2f | %7.2f |" %(c_cz, c_lty, c_rt))
    except ZeroDivisionError as e:
        print( s1, e)

print("--------------------------------------------------------------------------------------------------------------------------------------------------------------------")
avg_cz  = avg_cz  / len(benchresult)
avg_lty = avg_lty / len(benchresult)
avg_rt  = avg_rt  / len(benchresult)
print("==> avg:                                                                                                                            || %7.2f | %7.2f | %7.2f |" %(avg_cz, avg_lty, avg_rt))
print()
