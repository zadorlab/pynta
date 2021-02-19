import sqlite3
import os
from pathlib import Path
from ase.io import read
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo


class DataBase():

    def main(self, db_file, facetpath):
        sql_create_minima_table = ''' CREATE TABLE IF NOT EXISTS minima (
                                        id integer PRIMARY KEY,
                                        chemical_symbol text NOT NULL,
                                        prefix text,
                                        file_name text NOT NULL,
                                        total_energy real,
                                        zpe_energy real
                                    ); '''
        # create a database connection
        conn = self.create_connection(db_file)

        prepare = PrepareDataToDB(facetpath)

        minima_details = prepare.prepare_minima()
        minima_vib_details = prepare.add_zpe_energies()

        if conn is not None:
            self.create_table(conn, sql_create_minima_table)
        else:
            print('Error: cannot create the database connection')

        with conn:
            for _, minima in minima_details.items():
                self.create_minima(conn, minima)

            for key, value in minima_vib_details.items():
                prefix, chemical_symbol = key.split('_')
                vib_entry = (value, prefix, chemical_symbol)
                self.add_minima_vib(conn, vib_entry)
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
        sql = ''' INSERT INTO minima(chemical_symbol,prefix,file_name,total_energy)
              VALUES(?,?,?,?) '''
        cursor = conn.cursor()
        cursor.execute(sql, minima)
        conn.commit()
        # return cursor.lastrowid

    def add_minima_vib(self, conn, minima_vib):
        sql = ''' UPDATE minima
                SET zpe_energy = ?
                WHERE prefix = ? AND chemical_symbol = ?
                '''
        cursor = conn.cursor()
        cursor.execute(sql, minima_vib)
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
            prefix = os.path.basename(traj).split('.')[0]
            potential_energy = atoms.get_potential_energy()

            entry = (species, prefix, fname, potential_energy)

            details[i] = entry

        return details

    def add_zpe_energies(self):
        final_dict = {}
        path_to_vib_species = os.path.join(
            self.current_dir, self.facetpath, 'minima_vib')
        keyword = '**/*traj'
        trajs = Path(path_to_vib_species).glob(keyword)
        for traj in trajs:
            traj = str(traj)
            if 'vib.' not in traj:
                zpe_energy = self.get_zpe_energy(traj)
                final_dict.update(zpe_energy)
        return final_dict

    def get_zpe_energy(self, path_to_vib_species):
        zpe_energy_dict = {}
        atoms = read(path_to_vib_species)
        vib_path = os.path.dirname(path_to_vib_species)
        os.chdir(vib_path)

        prefix, species = os.path.basename(
            path_to_vib_species).split('.')[0].split('_')
        indices = [atom.index for atom in atoms if atom.position[2]
                   > atoms.cell[2, 2]/2.]
        vib = Vibrations(atoms, indices=indices)
        vib_energies = vib.get_energies()
        key = prefix + '_' + species
        try:
            thermo = IdealGasThermo(vib_energies, atoms)
            zpe_energy = thermo.get_ZPE_correction()
            zpe_energy_dict[key] = zpe_energy

        except ValueError:
            pass
        return zpe_energy_dict
