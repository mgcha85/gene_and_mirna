# import os
# import pandas as pd
#
#
# dirname = "D:/Users/mgcha/Dropbox/Studymaterials/Graduate/UCF/Courswork/CDA5106/project/MachineProblem2/sim_report"
# with open(os.path.join(dirname, 'myout0.txt')) as f:
#     lines = f.read().split('\n')
#
# report = pd.DataFrame()
# get = False
# m = benchmark = None
# for line in lines:
#     if './sim' in line:
#         m, benchmark = line.split(' ')[3:]
#         report.loc[m, benchmark.split('.')[0]] = None
#         get = True
#     if get and 'misprediction rate:' in line:
#         mrate = line.split('		')[1]
#         report.loc[m, benchmark.split('.')[0]] = float(mrate[:-1])
#         get = False
# report.to_excel(os.path.join(dirname, 'bimodal.xlsx'))

# import os
# import pandas as pd
#
#
# dirname = "D:/Users/mgcha/Dropbox/Studymaterials/Graduate/UCF/Courswork/CDA5106/project/MachineProblem2/sim_report"
# with open(os.path.join(dirname, 'myout1.txt')) as f:
#     lines = f.read().split('\n')
#
# report = {}
# get = False
# m = n = benchmark = None
# for line in lines:
#     if './sim' in line:
#         m, n, benchmark = line.split(' ')[3:]
#         bench = benchmark.split('.')[0]
#         report[bench] = report.get(bench, pd.DataFrame())
#         get = True
#     if get and 'misprediction rate:' in line:
#         mrate = line.split('		')[1]
#         report[bench].loc[m, n] = mrate[:-1]
#         get = False
#
# writer = pd.ExcelWriter(os.path.join(dirname, 'gshare.xlsx'), engine='xlsxwriter')
# for bench, df in report.items():
#     df.to_excel(writer, sheet_name=bench)
# writer.save()
# writer.close()

import os
import pandas as pd
import matplotlib.pyplot as plt

dirname = "D:/Users/mgcha/Dropbox/Studymaterials/Graduate/UCF/Courswork/CDA5106/project/MachineProblem2/sim_report"
fpath = os.path.join(dirname, "gshare.xlsx")
xls = pd.ExcelFile(fpath)
for sname in xls.sheet_names:
    df = xls.parse(sname, index_col=0)
    for idx in df.index:
        df.loc[idx].plot()
        # plt.show()
    plt.legend()
    plt.xlabel("n")
    plt.ylabel("miss rate (%)")
    plt.grid()
    # plt.show()
    plt.savefig(os.path.join(dirname, "figures/gshare_{}.png".format(sname)))
    plt.close()
