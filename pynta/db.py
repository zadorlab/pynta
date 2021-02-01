import sqlite3


class DataBase():
    def main(self, db_file):
        sql_create_projects_table = ''' CREATE TABLE IF NOT EXISTS projects (
                                        id integer PRIMARY KEY,
                                        name text NOT NULL,
                                        begin_date text,
                                        end_date text
                                    ); '''
        # create a database connection
        conn = self.create_connection(db_file)

        if conn is not None:
            self.create_table(conn, sql_create_projects_table)
        else:
            print('Error: cannot create the database connection')

        with conn:
            project = ('Testing sqlite3 with python',
                       '2021-01-01', '2021-02-02')
            project_update = ('It is working', '2021-02-02', '2021-01-01', 9)
            project_delete = (5,)
            # project_id = self.create_project(conn, project)
            # print(project_id)

            self.update_project(conn, project_update)
            self.delete_project(conn, project_delete)
            self.select_all_projects(conn)

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

    def create_project(self, conn, project):
        sql = ''' INSERT INTO projects(name,begin_date,end_date)
              VALUES(?,?,?) '''
        cursor = conn.cursor()
        cursor.execute(sql, project)
        conn.commit()
        return cursor.lastrowid

    def update_project(self, conn, project):
        sql = ''' UPDATE projects
                SET name = ?,
                    begin_date = ?,
                    end_date = ?
                WHERE id = ?
                '''
        cursor = conn.cursor()
        cursor.execute(sql, project)
        conn.commit()

    def select_all_projects(self, conn):
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM projects')

        rows = cursor.fetchall()

        for row in rows:
            print(row)

    def delete_project(self, conn, query):
        sql = '''
                DELETE FROM projects
                WHERE id = ?
                '''
        cursor = conn.cursor()
        cursor.execute(sql, query)
        conn.commit()
