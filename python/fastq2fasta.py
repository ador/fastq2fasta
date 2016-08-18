# Functions to transform fastq files to fasta files
# by cutting low quality reads (marked with # in the quality line)
# in our sample files with BASE=33, see  http://drive5.com/usearch/manual/quality_score.html

import sys
from datetime import datetime

MinSeqLen = 60
TrimBadEndings = True

def transformHeaderSimple(headerStr):
  return ">" + headerStr[1:]

def transformHeader(headerStr, n, origLen):
  pos = headerStr.rfind('=' + str(origLen))
  newHeader = ">" + headerStr[1:pos +1] + str(origLen - n)
  return newHeader

def cutSeq(seq, n):
  until = len(seq) - n
  return seq[0:until]

def convertFastqItem2FastaItem(inputLines, minSeqLen = MinSeqLen, trimBadEndings = TrimBadEndings):
  # there should be 4 input lines
  if (len(inputLines) != 4 or inputLines[0][0] != '@'):
    return "ERROR: WRONG INPUT FORMAT"
  header = inputLines[0]
  sequence = inputLines[1]
  qualityStr = inputLines[3]
  origLen= len(sequence)
  cnt = 0
  for c in reversed(qualityStr):
    if c != '#':
      break
    else:
      cnt += 1
  if (origLen - cnt < minSeqLen):
    return "TOO LOW SEQUENCE QUALITY : " + header
  if (cnt > 0 and trimBadEndings):
    return [transformHeader(header, cnt, origLen), cutSeq(sequence, cnt)]
  else:
    return [transformHeaderSimple(header), sequence]

def convertFastqFile2FastaFile( inFileName, outFileName, logFileName,
                                minSeqLen = MinSeqLen,
                                trimBadEndings = TrimBadEndings):
  infile = open(inFileName, 'r')
  outfile = open(outFileName, 'w')
  logfile = open(logFileName, 'w')
  # read input file in 4-line batches
  itemCounter = 0
  okCounter = 0
  linesToProcess = []
  for line in infile:
    linesToProcess.append(line.strip())
    if (len(linesToProcess) == 4):
      itemCounter += 1
      result = convertFastqItem2FastaItem(linesToProcess, minSeqLen, trimBadEndings)
      if (len(result) != 2):
        logfile.write(str(datetime.now()) + "  -  " + result + "\n")
      else:
        okCounter += 1
        for resLine in result:
          outfile.write(resLine + '\n')
      linesToProcess = []
  
  print("Items processed: " + str(itemCounter))
  infile.close()
  outfile.close()
  logfile.close()


def main(params):
  if len(params) < 3:
    print("Usage:")
    print("python fastq2fasta.py <input fastq file>  <output fasta file>  <log file>  [seq length threshold (default:60)] [trim low quality ends (default: True)]")
    return

  infile = params[0]
  outfile = params[1]
  logfile = params[2]
  if len(params) >= 5:
    trimEnd = bool(params[4])
    convertFastqFile2FastaFile(infile, outfile, logfile, minSeqLen, trimEnd)
  elif len(params) == 4:
    minSeqLen = int(params[3])
    print("Minimum sequence length is " + str(minSeqLen))
    convertFastqFile2FastaFile(infile, outfile, logfile, minSeqLen)
  else:
    print("Minimum sequence length is " + str(MinSeqLen))
    convertFastqFile2FastaFile(infile, outfile, logfile)
  # convertFastqFile2FastaFile('../sample_data/sample1.fastq',
  #      '../sample_data/sample1_result.fasta',
  #      '../sample_data/sample1_logs.txt')

if __name__ == "__main__":
  main(sys.argv[1:])
