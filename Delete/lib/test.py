import DOMysql

mysql = DOMysql.SQL("/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/Delete/config/config.txt")
insert_list = ["project_id","delivery_start_time","cloud_address","analysis_bool","delete_type","delivery_account","delivery_db_id","create_time"]
value_list = ['PM-XS01KF2024050238-02','2024-05-29 14-48-43', 'obs://annoroad-cloud-product/user/project/test/XS01KF2024050238/PM-XS01KF2024050238-02/ANNO_XS01KF2024050238_PM-XS01KF2024050238-02_2024-05-29/',0,'交付删除','others','101570','2024-05-29 14-48-43']

table = 'tb_deletion_info_test'
mysql.insert( table, insert_list, value_list)
