#!/usr/bin/env python3
import sys
import os
import argparse
import xlrd

'''
This software is available uner the MIT license:

Copyright (c) 2015 Sanjuro jogdeo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

'''

def main():
    args=processArgs()
    
    sys.stderr.write('\nProcessing {}: \n\n'.format(args.xlsfile))
    
    #get the data from the excel spreadsheet
    tableData, ignoredList = loadSpreadsheet(args.xlsfile)
    
    #validate that the required field attributes are specified and that the foreign key constraints are valid
    errors = validateTables(tableData)
    if len(errors) > 0 :
        for tbl, terr in errors.items() :
            if len(terr) != 0 :
                sys.stderr.write('Errors in the {} table: \n'.format(tbl))
                sys.stderr.write("\n".join([e.__str__() for e in terr]))
    
        sys.stderr.write("Exiting...")
        sys.exit()
    
    #print which columns in the spreadsheet have been ignored    
    for tablename,ilist in ignoredList.items() :
        sys.stderr.write('Ignored these columns in the {} table: {}\n'.format(tablename,','.join(ilist)))
    sys.stderr.write('\n\n')    
    
    #determine the order in which the tables should be created.  It matters because of the foreign keys
    allForeignKeys = {}
    tableOrder = []
    
    #get all foreign keys per table
    tableList = tableData.keys()
    for tablename,table in tableData.items() :
        for field in table.getFields() :
            if field.foreign_key :
                ftbl,ffield = field.foreign_key.split(':')
                if tablename in allForeignKeys :
                    allForeignKeys[tablename].append(ftbl)
                else :
                    allForeignKeys[tablename] = [ftbl]
    
    #get a list of tables without foreign keys                    
    for tname in tableList :
        if not tname in allForeignKeys :
            tableOrder.append(tname)
    
    #get the table order
    tableOrder = makeTableOrder(tableOrder,allForeignKeys,tableList,len(tableList)+1,0)
    if len(tableOrder) != len(tableList) :
        sys.stderr.write('Couldn\'t process the foreign keys. You might have specific \
                         cross-referenced keys, or this program might just not be up to the design. \n \
                         Exiiting.'
        )
    
    
    
    
    #make the sql
    sys.stderr.write('Writing the SQL statement\n\n')
    sql = makeMysql(tableData,tableOrder,args.dbname)
    sys.stdout.write(sql)


##################### - END MAIN - #####################


#generate the SQL for database creation based on the table data read from the spreadsheet
def makeMysql (tables,tOrder,dbname) :
    
    sql = 'CREATE DATABASE IF NOT EXISTS {};\n'.format(dbname)
    
    for tablename in tOrder :
        table = tables[tablename]
        pKey = ''
        sql += 'CREATE TABLE IF NOT EXISTS `{}`.`{}`(\n'.format(dbname,tablename)
        for field in table.getFields() :    
            sql += field.to_mysql()+',\n'
            if field.primary_key :
                pKey = field.fname
        sql += 'PRIMARY KEY ({})'.format(pKey)
        
        fKeys = table.getForeignKeys()
        if len(fKeys) == 0 :
            sql += '\n'
        else :
            sql += ',\n'
            
        fKeyLines = []
        for fKey in fKeys :
            dfield,fKeyText = fKey
            ftbl,ffield = fKeyText.split(':')
            fKeyLines.append('FOREIGN KEY ({}) REFERENCES {}({})'.format(dfield,ftbl,ffield))
        sql += ",\n".join(fKeyLines)
            
            
            
        sql += '\n);\n\n'
    
    
    return sql

#mysql will error if a foreign key is designated before it's referenced table is.  This
#routine recursively orders the tables based on whether their foreign key references have
#been defined or not.
def makeTableOrder(tOrder,fKeys,tList,maxLoops,loopCount) :
    if len(tOrder) >= len(tList) :
        sys.stderr.write('Tables have been fully ordered.\n\n')
        return tOrder
    elif maxLoops <= loopCount :
        sys.stderr.write('Hit the maximum number of table ordering loops.  Ordering probably didn\'t work.\n\n')
        return tOrder
    else :
        newKeys = {}
        for dtbl in fKeys :
            coveredTables = 0
            for ftbl in fKeys[dtbl] :
                if ftbl in tOrder or ftbl == dtbl :
                    coveredTables += 1
            if coveredTables == len(fKeys[dtbl]) :
                tOrder.append(dtbl)
            else :
                newKeys[dtbl] = fKeys[dtbl]
        
        loopCount += 1
        return makeTableOrder(tOrder,newKeys,tList,maxLoops,loopCount)
        
        
#load the table definitions from Excel
def loadSpreadsheet (xlsfilename) :
    xlObj = xlrd.open_workbook(xlsfilename)
    tables = {}
    
    aliases = Field.get_aliases()
    
    tables = {}
    ignored = {}
    for sheet in xlObj.sheets() :
        
        #check if it's a tab that we should be processing
        try :
            if sheet.cell_value(0,0) != 'field name' :
                continue
        except IndexError as e :
            continue
        
        tableName = sheet.name

        #first get the column number associated with each field attribute.  This allows for each tab to
        #list different sets of columns and to have them in a different order
        fColumns = {}
        for i in range(0,len(sheet.row_values(0))) :
            fieldOnSheet = sheet.cell_value(0,i)
            if Field.getFieldFromAlias(fieldOnSheet) :
                fColumns[Field.getFieldFromAlias(fieldOnSheet)] = i
            else :
                if tableName in ignored :    
                    ignored[tableName].append(fieldOnSheet)
                else :
                    ignored[tableName] = [fieldOnSheet]
                

        #now create the table definition        
        t = Table()

        #loop through rows and save field attribs to field objects
        for r in range(1,sheet.nrows) :
            
            #Don't continue processing if there is a blank row
            row_values = sheet.row_values(r)            
            for val in row_values :
                if val != '' :
                    break
            else :
                break
                
            f=Field()
            for attrib in fColumns :
                value = sheet.cell_value(r,fColumns[attrib])        
                f.set_fieldAttribs(**{attrib:value})
            t.addField(f)
      
        tables[tableName] = t
        
    return tables,ignored
    
#makes sure the table definitions pass certain checks
def validateTables (tables) :
    foreignKeys = []
    vErrors = {}
    for tablename,table in tables.items() :
        
        #first run field validation
        elist = table.validate()
        if len(elist) > 0 : 
            vErrors[tablename] = table.validate()
        
        #now get foreign keys for use in the next step
        for field in table.getFields() :
            if field.foreign_key :
                foreignKeys.append((tablename,field.foreign_key))

    for tablename,fkstr in foreignKeys :
        
        try :
            fkList = fkstr.split(':')
            if len(fkList) != 2 :
                raise TableException('A foreign key specified as {} in the {} table is not valid'.format(fkstr,tablename))
            ftbl,ffield = fkList
            test = ffield in tables[ftbl].getFieldNames()
        except Exception as e :
            
            te = TableException('A foreign key specified as {} in the {} table is not valid'.format(fkstr,tablename))
            if tablename in vErrors.keys() :
                vErrors[tablename].append(te)
            else :
                vErrors[tablename] = [te]
            
    return vErrors
            
        
    


class Field (object):
    
    
    #set class variables
    
    requiredAttribs = ['fname', 'ftype']
    
    #aliases allow for readability and flexibility of columnn names in the excel spreadsheet
    fieldAliases = {
        'fname':['field name','field_name'],
        'ftype':['field type','field_type'],
        'unique': ['unique'],
        'autoincrement':['auto increment','autoincrement'],
        'primary_key':['primary key'],
        'foreign_key':['foreign key'],
        'not_null':['not null', 'not_null']
    }
    fieldFromAlias = {}
    for f, aList in fieldAliases.items() :
        for a in aList :
            fieldFromAlias[a]=f
    
    ftypeAliases = {
        'INT':['int','integer'],
        'DATE':['date'],
        'TEXT':['text'],
        'FLOAT':['float','decimal'],
        
    }
    
    ftypeFromAlias = {}
    for f, aList in ftypeAliases.items() :
        for a in aList :
            ftypeFromAlias[a]=f
    
    
    #set init and other methods        
    def __init__(self,**kwargs) :   
        self.fname = False
        self.ftype = False
        self.uniqueFlag = False
        self.autoincrement = False
        self.primary_key = False
        self.foreign_key = False
        self.not_null = False
        
        self.set_fieldAttribs(**kwargs)
        
    
    @staticmethod
    def get_aliases() :
        return Field.fieldAliases
    
    @staticmethod
    def get_fieldList() :
        return Field.fieldAliases.keys()
    
    @staticmethod
    def getFieldFromAlias(alias) :
        try :
            return Field.fieldFromAlias[alias]
        except KeyError as ke :
            return False
        except Exception as e :
            sys.stderr.write(e)
            sys.stderr.write('Catastrophic error. Exiting.\n')
            sys.exit()
    
    @staticmethod
    def getftypeFromAlias(alias) :
        try :
            return Field.ftypeFromAlias[alias]
        except KeyError as ke :
            return False
        except Exception as e :
            sys.stderr.write(e)
            sys.stderr.write('Catastrophic error. Exiting.\n')
            sys.exit()
    
    
    def set_fieldAttribs(self,**kwargs) :
        setErrors = []
        for arg in kwargs :
            if arg.lower() in Field.fieldAliases :
                if kwargs[arg] != '' :
                    setattr(self,arg,kwargs[arg])
                else :
                    arg == False
                
            else :
                setErrors.append('the column name \'{}\' is not understood by this script \
                                 as a valid database field attribute'.format(arg))
        
        return setErrors
        
    def validate (self) :
        
        exceptionText = ''
        
        #check if the required fields are present
        validationErrors = []
        priorFtypeError = False
        for rAttrib in Field.requiredAttribs :
            atval=getattr(self,rAttrib)
            if getattr(self,rAttrib) == False :
                validationErrors.append(rAttrib)
                if rAttrib == 'ftype' :
                    priorFtypeError = True
        if len(validationErrors) > 0 :
            exceptionText = "A field is missing the following attribs: {}\n".format(",".join(validationErrors))
        
        #check if ftype is a recognized value
        if not Field.getftypeFromAlias(self.ftype) and not priorFtypeError :
            exceptionText += "A field type of {} is not recognized as a valid type\n".format(self.ftype)
        
        if exceptionText != '' :
            raise TableException(exceptionText)
        else :
            return True
        
    def to_mysql (self) :
        
        defList = [self.fname]
        defList.append(Field.getftypeFromAlias(self.ftype))
        if self.uniqueFlag :
            defList.append('UNIQUE')
        if self.autoincrement :
            defList.append('AUTO_INCREMENT')
        if self.not_null :
            defList.append('NOT NULL')
        
        return ' '.join(defList)
            
    
class TableException(Exception) :
    pass

                                                                 
class Table :
    
    def __init__(self) :
        self.fields = []
        
    def addField (self,fieldObj) :
        
        if isinstance(fieldObj,Field) :
            self.fields.append(fieldObj)
    
    def getFields(self) :
        return self.fields
    
    def getFieldNames(self) :
        names = []
        for field in self.getFields() :
            names.append(field.fname)
        return names
    
    def getForeignKeys(self) :
        fKeys = []
        for field in self.getFields() :
            if field.foreign_key :
                fKeys.append((field.fname,field.foreign_key))
        return fKeys
            
    
    def validate (self) :
        tValidationErrors = []
        
        for field in self.fields :
            try :
                field.validate()
            except TableException as fe :
                tValidationErrors.append(fe)
        
        return tValidationErrors 


def processArgs():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('\nError: %s\n\n' % message)
            self.print_help()
            sys.exit(2)
    
    class CheckFile(argparse.Action) :
        def __call__(self,parser,namespace,value,option_string) :
            if os.path.isfile(value)==False :
                parser.error("The excel_file argument flag needs a valid filename")
                parser.print_help()
            else :
                setattr(namespace,self.dest,value)
            
    
    
    
    #argParser = MyParser(usage=("%s (sourceDir & filter) | filterFile" % (os.path.basename(sys.argv[0]))))
    argParser = MyParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                         description='''
Converts Excel spreadsheet to SQL for database creation.

Individual tabs represent tables. The first row of each worksheet should
contain field attribute headers.  Any worksheet that doesn't contain
the text 'field_name' in cell A1 will be ignored.  A blank row will
end field imports.  Foreign keys should be designated by the referenced table
and reference field with a colon in between (i.e. tbl:field).

A very limited number of field types are currently supported (int, float,
date, and text).  


'''
                         )
    
    argParser.add_argument('dbname', metavar="dbname", help="The name of the database")
    argParser.add_argument('xlsfile', metavar="excel_file", action=CheckFile, help="The Excel schema file")
    
    ap=argParser.parse_args()
    return ap



    
#This is required because by default this is a module.  Running this makes it execute main as if it is a script
if __name__ == '__main__':
    main()
    
