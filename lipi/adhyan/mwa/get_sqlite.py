import sqlite3
import pandas as pd

# creating file path
dbfile = '/media/rohit/MWA/cigale/data.db'
# Create a SQL connection to our SQLite database
con = sqlite3.connect(dbfile)

# creating cursor
cur = con.cursor()

# reading all table names
table_list = [a for a in cur.execute("SELECT name FROM sqlite_master WHERE type = 'table'")]
# here is you table list
print(table_list)
bc03=pd.read_sql_query('SELECT * FROM bc03', con)
filters=pd.read_sql_query('SELECT * FROM filters', con)

def importdb(db):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("SELECT name FROM sqlite_master WHERE type='table';")
    for table in c.fetchall():
        yield list(c.execute('SELECT * from ?;', (table[0],)))

dd=importdb(dbfile)

# Be sure to close the connection
con.close()

