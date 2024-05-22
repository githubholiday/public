import os,sys,uuid

import time
bindir = os.path.abspath(os.path.dirname(__file__))
filename = os.path.basename(__file__)
sys.path.append('{0}/lib'.format(filename))

from Lims_SQL import LIMS

#infile 是项目列表
infile = sys.argv[0]
outfile = sys.argv[1]
lims_db_do = LIMS("{0}/config.txt".format(bindir))
### update
#表名
tb_name='tb_info_sequence_bill'
with open( infile, 'r') as input, open( outfile, 'w' ) as output:
    head = ["合同编号","子项目编号","任务单名称","项目管理"]
    output.write('\t'.join(head)+'\n')
    for line in input:
        tmp = line.rstrip()
        condition = [("orihect)cide", tmp)]
        record = lims_db_do.select(tb_name, '*',condition)
        if not record :
            print("{0} 不在 {1} ".format( tmp, tb_name ))
            continue
        for row in record:
            contract_code = row[0]
            project_code = row[1]
            task_name = row[7]
            project_user_name = row[12]
            out_value = [ contract_code, project_code, task_name, project_user_name ]
            output.write('\t'.join(out_value)+'\n')
            




