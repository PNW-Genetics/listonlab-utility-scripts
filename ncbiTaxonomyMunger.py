#!/bin/env python2.7
import sys
import os
import argparse
from collections import defaultdict
import pickle

'''
License:

This is free and unencumbered software released into the public domain.  All authors are or
were bona fide officers or employees of the U.S. Government at the time the software was
developed and that the work is a “work of the U.S. Government” (prepared by an officer or
employee of the U.S. Government as a part of official duties), and, therefore, is not subject
to U.S. copyright as stated in Patent, Trademark and Copyright Laws 17 U.S.C. §105.
 
The software is provided “as is”, without warranty of any kind, express or implied, including
but not limited to the warranties of merchantability, fitness for a  particular purpose and
non-infringement. In no event shall the authors be liable for any claim, damages or other
liability, whether in an action of contract, tory or otherwise, arising from, out of or in
connection with the software or the use of other dealings in the software.
'''

def main():

    args=processArgs()

    rankOrder = ['superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum',
             'superclass', 'class', 'subclass', 'infraclass',
             'superorder', 'order', 'suborder', 'parvorder', 'infraorder',
             'superfamily', 'family', 'subfamily', 'tribe', 'subtribe', 'genus', 'subgenus',
             'species group', 'species subgroup', 'species', 'subspecies', 'varietas', 'forma']
    revRanks = list(rankOrder).reverse()
    majorRanks = ['kingdom', 'phylum','class', 'order', 'family', 'genus', 'species']
    revMajor = list(majorRanks).reverse()
    
    
    #load db either from a pickle or directly from db files
    if args.p != None :
        sys.stderr.write("Reading pickle...\n")
        with open(args.p,'rb') as lh :
            taxdb = pickle.loads(lh.read())
        
    else :
        #read taxonomy db files
        nodefile = os.path.join(args.t,'nodes.dmp')
        namefile = os.path.join(args.t,'names.dmp')
        
        #put names into dictionary keyed by taxid
        nameDict = readNames(namefile)
        sys.stderr.write("Names in names.dmp: %s\n" % (len(nameDict)))
        
        nodeDict,fileRanks = readNodes(nodefile)
        sys.stderr.write("Nodes in nodes.dmp: %s\n" % (len(nodeDict)))
        
        taxdb = builddb(nodeDict,nameDict,rankOrder)
        
        
        #pickle taxdb
        if args.s != None :
            with open (args.s,'wb') as ph:
                pickle.dump(taxdb,ph)
        
    
    #build the list of ranks to be output
    outRanks = ['kingdom']
    if args.r != None :
        if args.r == 'all' :
            outRanks = list(revRanks)
        elif args.r == 'major' :
            outRanks = list(revMajor)
        else :
            outRanks = args.r.split(',')
    
    
    
    #set the input and output configuration based on -i and -c flags
    if args.out == '-' :
        outfh = sys.stdout
    else :
        outfh = open(args.out,'w')
        
    if args.taxids == '-' :
        infh = sys.stdin
    else :
        infh = open(args.taxids,'r')
    
    
    #add sciname to the output ranks as standard practice
    expOutRanks = list(outRanks)
    expOutRanks.insert(0,'sciName')
    
    #need to tell user what columns are being added since we're not printing them yet.
    sys.stderr.write("These are columns that will be output:\n")
    sys.stderr.write("taxid\t"+"\t".join(expOutRanks)+"\n\n\n")
    
    
    
    skipLines = args.k or 0
    lineCounter = 1
    #loop through input lines/taxids and output taxonomy info
    discards = []
    for line in infh :
        
        if lineCounter <= skipLines :
            lineCounter += 1
            continue
        
            
        line = line.strip()
        
        #set default taxid and prefix, assuming only a signle column of taxids is entered
        taxid = getTaxID(line)
        prefix = "%s\t" % (line)
        
        #now handle situation where input is multi-column
        if args.c != None :
            prefix = line+"\t"
            lineParts = line.split("\t")
            try :
                taxid = getTaxID(lineParts[args.c-1])
            except IndexError:
                print "caught ",line
                discards.append(line)
                
            
        outfh.write(prefix)
        
        rankVals = list()
        for orank in expOutRanks :
            try :
                rankVals.append(taxdb[taxid][orank])
            except KeyError as e :
                rankVals.append('NA')
        outfh.write("\t".join(rankVals)+"\n")
        
    outfh.close()
    infh.close()
    
    if len(discards) > 0 :
        if len(discards) < 50 :
            sys.stderr.write("The following lines had fewer columns than he -c flag\n")
            sys.stderr.write("\n".join(discards)+"\n\n")
        else :
            sys.stderr.write("There were %s lines with fewer columns than he -c flag.  Here are the first 50:\n" % (len(discards)))
            sys.stderr.write("\n".join(discards[:50])+"\n\n")
    

def getTaxID (taxField) :

    taxids = taxField.split(';')
    try :
        rTaxid = int(taxids[0])
        return str(rTaxid)
    
    except :
        return "NA"
    

#builds the dict of dict database struction that will be used for looking up taxonomy ranks
def builddb (nodes,names,ranks) :
    sys.stderr.write("Starting db build...\n\n")
    indexes = list(ranks).insert(0,'sciName')
    taxDict = defaultdict(dict)
    
    stopper = 1
    for taxid in nodes :
        nodelist = []
        taxInfo = getParentInfo(taxid,nodes,ranks,nodelist)
        try :
            #this is in a try block in case the taxid wasn't captured for some reason.
            #Might happen with old pickles.  Prefer catastrophic failure for now.
            taxDict[taxid]['sciName']=names[taxid]
        except Error as e:
            sys.stderr.write("Error, exiting: %s\n" % (e))
            sys.exit()
            
        #for each rank returned, save into the dict using the rank itself (e.g. genus) as key,
        #and the rank value (e.g. Homo) as the value.
        for pid,rank in taxInfo :
            try :
                pname = names[pid]
                taxDict[taxid][rank]=names[pid]
                
            except Error as e:
                sys.stderr.write("Error, exiting: %s\n" % (e))
                sys.exit()
        
        #output progress
        stopper += 1
        if stopper % 1000 == 0 :
            sys.stderr.write("%s processed\r" % (stopper))
        #if stopper >= 3000 : break
        
    return taxDict
    
    
#recursive walk up the parent tree
def getParentInfo (taxid,nodes,ranks,nodelist) :
    
    pid,rank = nodes[taxid][0:2]
    
    if nodes[taxid][0] == '1' :
        if rank in ranks : 
            nodelist.append((taxid,nodes[taxid][1]))
        return nodelist
    else :
        if rank in ranks :  
            nodelist.append((taxid,rank))
          
        getParentInfo(nodes[taxid][0],nodes,ranks,nodelist)
        return nodelist




def readNodes (nodefile) :
    nodeDict = dict()
    uniqRanks = dict()
    
    with open(nodefile,'r') as fh :
        for line in fh :
            taxid,parentid,rank =  line.split('|')[0:3]
            taxid = taxid.strip(' \t')
            
            rank = rank.strip(' \t')
            uniqRanks[rank] = 1
            
            parentid = parentid.strip(' \t')
            
            nodeDict[taxid] = (parentid,rank)
            
    return nodeDict,uniqRanks
    
    
def readNames(namefile) :
    nameReturn = dict()
    with open(namefile,'r') as fh :
        for line in fh :
            
            #don't need to load all the extra names since there is always at least one
            #scientific name per taxid
            if 'scientific name' in line :
                taxid,taxname =  line.split('|')[0:2]
                taxid = taxid.strip(' \t')
                taxname = taxname.strip(' \t')
                nameReturn[taxid] = taxname
    return nameReturn
                
    

def processArgs(helpFlag=False):
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('\nError: %s\n\n' % message)
            self.print_help()
            sys.exit(2)
    
    class CheckDBFiles(argparse.Action) :
        def __call__(self,parser,namespace,value,option_string) :
            nodefile = os.path.join(value,'nodes.dmp')
            namefile = os.path.join(value,'names.dmp')
            exitFlag = 0
            if os.path.isfile(nodefile) == False :
                sys.stderr.write("Couldn't find nodes.dmp file.\n")
                exitFlag = 1
            if os.path.isfile(namefile) == False :
                sys.stderr.write("Couldn't find names.dmp file.\n")
                exitFlag = 1
            if exitFlag == 1 :
                parser.error("Couldn't find the taxonomy names.dmp and/or nodes.dmp files")
                
            else :
                setattr(namespace,self.dest,value)
    class CheckFileExists(argparse.Action) :
        def __call__(self,parser,namespace,value,option_string) :
            if value != '-' and os.path.isfile(value) == False :
                parser.error("Couldn't find the file specified by %s: %s\n" % (option_string,value))
            else :
                setattr(namespace,self.dest,value)

    
    
    
    #argParser = MyParser(usage=("%s (sourceDir & filter) | filterFile" % (os.path.basename(sys.argv[0]))))
    argParser = MyParser(description="""
    Accepts a list of NCBI taxids and outputs taxonomic information from a local
    copy of the NCBI taxonomy database.  Taxids can be read from STDIN, or a file.
    If a tab delimited file provided, taxonomic information will be appended to
    the end of each row.  In this case, taxids must be provided in one of the
    columns of the tab delimited file.  Blast output for taxid often contains a
    list semicolon separated taxids.  The script uses the first taxid in the list
    to search for the rest of the taxonomy information.  
    
    
    User must provide either a pickled database file or a path to the NCBI taxonomy 
    database directory.  This directory must contain names.dmp and nodes.dmp files.
    Surprisingly, pickles take about twice as long to load, but they can be passed
    around and are a fixed record of the database used.  As you might expect though,
    old pickles don't taste very good.
    
    
    User can specify which taxonomic ranks should be output using commma separated
    values in the -r flag, for instance, by specfiying 'family,kingdom'.  Specifying
    'major' will output the seven major taxonomic ranks.  Specifying 'all' will yield
    a column for each of the 28 named ranks, though many of these are likely to be
    null values.  The raw taxonomy database contains many levels designated as
    'no rank'.  These are never output by this script.
    
                         """,
                         formatter_class=argparse.RawTextHelpFormatter)
    
    argParser.add_argument('taxids', metavar='taxids', action=CheckFileExists, help='A file containing taxids.  Alternately, use \'-\' to take values from STDIN.')
    argParser.add_argument('-c', metavar='int (column)',  type=int, help='Used when taxids are a column ' \
                           'within a tab separated file. The column number should be indicated with the -c flag.')
    
    argParser.add_argument('out', metavar="outname", help="Output filename, or use '-' for STDOUT.")
    
    dbinputs = argParser.add_mutually_exclusive_group(required=True)
    dbinputs.add_argument('-p', metavar='file', help='A pickled taxonomy database.')
    dbinputs.add_argument('-t', metavar='directory', action=CheckDBFiles,  help='Taxonomy db directory')
    
    argParser.add_argument('-s', metavar="filename", help="Serialize the taxonomy database structure and save in the file specified.")
    argParser.add_argument('-k', metavar="int", type=int, help="Number of rows at top to skip.")
    argParser.add_argument('-r', metavar="ranks", help='Comma separated list of ranks to output. ' \
                           'Use \'all\' for all ranks. Default=scientific name, kingdom.  Sci name ' \
                           'will be output regardless of rank choices.')
    
    
    if helpFlag==True :
        argParser.print_help()
        return
        
    else :
        ap=argParser.parse_args()
        return ap




#This is required because by default this is a module.  Running this makes it execute main as if it is a script
if __name__ == '__main__':
    main()
    