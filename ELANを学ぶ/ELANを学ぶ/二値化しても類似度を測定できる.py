#! C:\Program Files (x86)\Microsoft Visual Studio\Shared\Python36_64\
# -*- coding: utf-8 -*-
import os
import sys
import glob
import subprocess
import datetime # datetimeモジュールのインポート
# Python 2.7, 3.4, 3.6 が複数の場所にある場合、'py -3.6 xxx' のように指定して実行する
# 例：　py -3.6 -m pip install --upgrade pip
CURRENT_PATH = os.path.dirname(__file__)
EDITOR_PATH = '"C:\\Program Files (x86)\\sakura\\sakura.exe"'
# pip install mysqlclient
import MySQLdb as my
# DB connect test
con = my.connect(user='root', password='', database='atted1388data')#, use_unicode=True) #, charset="utf8")
cursor = con.cursor()
statement = "select * from expression_data_part2 where AGI = 'At2g26470';"
cursor.execute(statement)
records = cursor.fetchall()
con.close()
# record を確認 
for record in records:
    print(record)
