import os
import sqlite3
import pandas as pd
import numpy as np
from dateutil.rrule import rrule, MONTHLY, YEARLY
import datetime


class Database:
    def __init__(self, root, fname, con=False):
        self.root = root
        self.fname = fname
        if con is True:
            self.con = sqlite3.connect(os.path.join(root, 'database', fname + '.db'), check_same_thread=False)

    def load_db(self, fname):
        out_path = os.path.join(self.root, 'database/{}'.format(fname))
        return sqlite3.connect(out_path, check_same_thread=False)

    def get_dairly_data(self, code, date, interval='month'):
        if interval == 'month':
            con = sqlite3.connect(os.path.join(self.root, self.fname + '_{}-{}.db'.format(date.year, date.month)), check_same_thread=False)
        elif interval == 'year':
            con = sqlite3.connect(os.path.join(self.root, self.fname + '_{}.db'.format(date.year)), check_same_thread=False)
        else:
            return
        df = pd.read_sql_query("SELECT * FROM {tn} WHERE Date BETWEEN '{date}' AND '{date}'".format(tn=code, date=date), con, index_col='Date')
        df.index = pd.to_datetime(df.index)
        return df

    def get_data_by_time(self, code, date, time1, time2, interval='month'):
        if interval == 'month':
            con = sqlite3.connect(os.path.join(self.root, self.fname + '_{}-{}.db'.format(date.year, date.month)), check_same_thread=False)
        elif interval == 'year':
            con = sqlite3.connect(os.path.join(self.root, self.fname + '_{}.db'.format(date.year)), check_same_thread=False)
        else:
            return

        df = pd.read_sql_query("SELECT * FROM {tn} WHERE Date BETWEEN '{date} {time1}' AND '{date} "
                               "{time2}'".format(tn=code, date=date, time1=time1, time2=time2), con, index_col='Date')
        df.index = pd.to_datetime(df.index)
        return df

    def get_range(self, date1, date2, interval='month'):
        if interval == 'month':
            dates = [dt.date() for dt in rrule(MONTHLY, dtstart=date1, until=date2)]
            if date2.day < date1.day:
                dates.append(date2)
            return ['{:04d}-{:02d}'.format(x.year, x.month) for x in dates]
        elif interval == 'year':
            dates = [dt.date() for dt in rrule(YEARLY, dtstart=date1, until=date2)]
            if date2.month < date1.month:
                dates.append(date2)
            return ['{:04d}'.format(x.year) for x in dates]

    def get_data_by_datetime(self, code, date1, date2, time1, time2, interval='month'):
        range = self.get_range(date1, date2, interval)
        dfs = []
        for date in range:
            con = sqlite3.connect(os.path.join(self.root, 'database', self.fname + '_{}.db'.format(date)), check_same_thread=False)
            if Database.checkTableExists(con, code) is False:
                return
            
            df = pd.read_sql_query("SELECT * FROM '{tn}' WHERE Date BETWEEN '{date1} {time1}' AND '{date2} {time2}'".format(tn=code,
                                    date1=date1, date2=date2, time1=time1, time2=time2), con, index_col='Date').astype(float)
            df.index = pd.to_datetime(df.index)
            dfs.append(df)
        return pd.concat(dfs)

    def get_data_by_date(self, code, date1, date2, interval='month'):
        range = self.get_range(date1, date2, interval)
        dfs = []
        for date in range:
            con = sqlite3.connect(os.path.join(self.root, 'database', self.fname + '_{}.db'.format(date)), check_same_thread=False)
            df = pd.read_sql_query("SELECT * FROM '{tn}' WHERE Date BETWEEN '{date1}' AND '{date2}'".format(tn=code,
                                date1=date1, date2=date2), con, index_col='Date')
            df.index = pd.to_datetime(df.index)
            dfs.append(df)
        return pd.concat(dfs)

    def get_data(self, code, interval='month'):
        today = datetime.date.today()

        if interval == 'month':
            con = sqlite3.connect(os.path.join(self.root, 'database', self.fname + '_{}-{}.db'.format(today.year, today.month)), check_same_thread=False)
        elif interval == 'year':
            con = sqlite3.connect(os.path.join(self.root, 'database', self.fname + '_{}.db'.format(today.year)), check_same_thread=False)
        else:
            return

        df = pd.read_sql_query("SELECT * FROM '{tn}'".format(tn=code), con, index_col='Date')
        df.index = pd.to_datetime(df.index)
        return df

    @staticmethod
    def get_years(dateList):
        years = []
        for date in dateList:
            date = '{:04d}'.format(date.year)
            if date not in years:
                years.append(date)
        return years

    @staticmethod
    def get_months(dateList):
        months = []
        for date in dateList:
            date = '{:04d}-{:02d}'.format(date.year, date.month)
            if date not in months:
                months.append(date)
        return months

    @staticmethod
    def get_dateList(root):
        flist = os.listdir(root)
        flist = [x for x in flist if 'intra_day_5min' in x]

        dateList = []
        for fname in flist:
            fpath = os.path.join(root, fname)
            con = sqlite3.connect(fpath, check_same_thread=False)
            df = pd.read_sql_query("SELECT * FROM '038500'", con, index_col='Date')
            df.index = pd.to_datetime(df.index)
            dateList.extend([*np.unique(df.index.date)])
        return dateList

    @staticmethod
    def get_sample(root):
        today = datetime.date.today()
        con = sqlite3.connect(os.path.join(root, 'database', 'intra_day_5min_{:04d}-{:02d}.db'.format(today.year, today.month)), check_same_thread=False)
        df = pd.read_sql_query("SELECT * FROM '038500'", con, index_col='Date')
        df.index = pd.to_datetime(df.index)
        return df

    @staticmethod
    def load_tableList(con):
        cursor = con.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()
        df_slist = []
        for table_name in tables:
            df_slist.append(table_name[0])
        return np.array(sorted(df_slist))

    @staticmethod
    def checkTableExists(con, tablename):
        dbcur = con.cursor()
        dbcur.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='{}'".format(tablename))
        if dbcur.fetchone()[0] == 1:
            dbcur.close()
            return True
        dbcur.close()
        return False

    @staticmethod
    def load_table_columns(con, tname):
        cursor = con.execute("SELECT * FROM '{}'".format(tname))
        return list(map(lambda x: x[0], cursor.description))
