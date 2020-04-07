# =============================================
# Example code: how to use the cbl.Table object
# =============================================

# import the CosmoBolognaLib modules 
import CosmoBolognaLib as cbl

# define the column names
names = ["c1", "c2", "c3", "c4", "c5"]

# define the index of the columns to be read
use_cols = [] 

# define the number of header lines
header_lines_to_skip = 2

# define the Table object, reading from file
table = cbl.Table("./", "data.dat", names, use_cols, header_lines_to_skip)

# print the column with name "c5"
print(table["c5"])

# insert a new column wit name c6
table.insert("c6", [-1, -2, -3])

# write the table on an ascii file
table.write("./", "data_out.dat")
