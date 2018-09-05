import csv
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument('--bim', help = 'file path of the target bim file to be translated')
parser.add_argument('--csv', help = 'file path of the target affymetrix csv file containing the translations')
parser.add_argument('-o', '--outputdir', help = 'Output directory for the new bim file')
args = parser.parse_args()

AffyDict={}
with open(args.csv, 'r') as affy:
  for line in affy:
    if ( line[0] != '#' ):
     affy.next()
    else:
      fields = line.split(',')
      if (fields[1] != '---'):
        AffyDict[fields[0]]=fields[1]

newbim = []
with open(args.bim, 'r') as bimfile:
  for line in bimfile:
    fields = line.split('\t')
    if ( fields[1] in AffyDict ):
      fields[1] = AffyDict[fields[1]]
    newbim.append('\t'.join(fields))
     
outbim = '/00AffyTranslated.bim'
outbim = args.outputdir + outbim
with open(outbim, 'wb') as out:
  for line in newbim:
    out.write(line + '\n')
      
      

      
