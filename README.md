Liston/Cronn Lab utility scripts
======

##bcsort_v5.pl##

A perl program to partition barcoded fasta or fastq files.  


##old_bcsort_versions##

Versions of bcsort saved separately, some are associated with various publications


##excelToSchema.py##

A python script for converting an Excel template to an executable SQL table creation script.
Used for the redevelopment of the wildstrawberry website and still a work in progress.  

Obvious next steps:
* support more field types
* support more types of foreign key constraints
* add validation to ensure that each table has a primary key
* tighten up whether boolean fields are considered True or not.  Currently entering 'False' into the spreadsheet could be interpreted as True just because it is not null
* add tabs that have been ignored to the ignore list
* add check to makes sure fields with a foreign key ref have an int type


