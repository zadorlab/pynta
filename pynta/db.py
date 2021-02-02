import sqlite3
import os
import io
import numpy as np
from pathlib import Path
from ase.io import read


class DataBase():

    def main(self, db_file, facetpath):
        sql_create_minima_table = ''' CREATE TABLE IF NOT EXISTS minima (
                                        id integer PRIMARY KEY,
                                        chemical_symbol text NOT NULL,
                                        file_name text NOT NULL,
                                        total_energy real
                                    ); '''
        # create a database connection
        conn = self.create_connection(db_file)

        minima_details = PrepareDataToDB(facetpath).prepare_minima()

        if conn is not None:
            self.create_table(conn, sql_create_minima_table)
        else:
            print('Error: cannot create the database connection')

        with conn:
            for _, minima in minima_details.items():
                self.create_minima(conn, minima)
            # minima = ('Testing sqlite3 with python',
            #         '2021-01-01', '2021-02-02')
            # minima_update = ('It is working', '2021-02-02', '2021-01-01', 9)
            # minima_delete = (5,)
            # print(minima_id)

            # self.update_minima(conn, minima_update)
            # self.delete_minima(conn, minima_delete)
            self.select_all_minima(conn)

    def create_connection(self, db_file):
        conn = None
        try:
            conn = sqlite3.connect(db_file)
        except sqlite3.Error as e:
            print(e)
        return conn

    def create_table(self, conn, create_table_sql):
        try:
            cursor = conn.cursor()
            cursor.execute(create_table_sql)
        except sqlite3.Error as e:
            print(e)

    def create_minima(self, conn, minima):
        sql = ''' INSERT INTO minima(chemical_symbol,file_name,total_energy)
              VALUES(?,?,?) '''
        cursor = conn.cursor()
        cursor.execute(sql, minima)
        conn.commit()
        # return cursor.lastrowid

    def update_minima(self, conn, minima):
        sql = ''' UPDATE minima
                SET name = ?,
                    begin_date = ?,
                    end_date = ?
                WHERE id = ?
                '''
        cursor = conn.cursor()
        cursor.execute(sql, minima)
        conn.commit()

    def select_all_minima(self, conn):
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM minima')

        rows = cursor.fetchall()

        for row in rows:
            print(row)

    def delete_minima(self, conn, query):
        sql = '''
                DELETE FROM minima
                WHERE id = ?
                '''
        cursor = conn.cursor()
        cursor.execute(sql, query)
        conn.commit()


class PrepareDataToDB():
    def __init__(self, facetpath):
        self.current_dir = os.getcwd()
        self.facetpath = facetpath
        self.path = os.path.join(self.current_dir, self.facetpath)

    def prepare_minima(self):
        details = {}
        minima_path = os.path.join(self.path, 'minima')
        keyword = '**/*traj'
        trajs = Path(minima_path).glob(keyword)
        for i, traj in enumerate(trajs):
            species = os.path.basename(os.path.dirname(traj))
            fname = os.path.join('./' + self.facetpath, 'mimima',
                                 os.path.basename(traj))
            atoms = read(traj)
            potential_energy = atoms.get_potential_energy()
            entry = (species, fname, potential_energy)
            details[i] = entry
        return details
