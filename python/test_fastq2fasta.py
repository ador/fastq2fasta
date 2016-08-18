from fastq2fasta import *

def test_read_fastq():
  inData = open('../sample_data/sample1.fastq', 'r')
  counter = 0
  for line in inData:
    counter += 1
  assert(200, counter)

def test_cut_seq():
  sequence = "123456"
  assert("1234" == cutSeq(sequence, 2))
  assert("123456" == cutSeq(sequence, 0))
  assert("" == cutSeq(sequence, 6))

def test_transform_header_simple():
  header = "@lorem ipsum dolor length=60"
  expHeader1 = ">lorem ipsum dolor length=60"
  assert(expHeader1 == transformHeaderSimple(header))

def test_transform_header():
  header = "@lorem ipsum dolor length=60"
  expHeader1 = ">lorem ipsum dolor length=49"
  assert(expHeader1 == transformHeader(header, 11, 60))
  expHeader2 = ">lorem ipsum dolor length=40"
  assert(expHeader2 == transformHeader(header, 20, 60))

def test_trim_bad_endings1():
  oneFastqItem = [ "@ERR445333.1.2 1 length=30",
        "AACCCCTAACCCCTAACCCTAACCCTAACC",
        "+ERR445333.1.2 1 length=30",
        "EFFFGF?EFGA???DC>:?AEEEAAAAAAA" ]
  resultLines = convertFastqItem2FastaItem(oneFastqItem, 10)
  assert(2, len(resultLines))
  assert('>ERR445333.1.2 1 length=30' == resultLines[0])
  assert(30, len(resultLines[1])) # nothing should be cut

def test_trim_bad_endings2():
  oneFastqItem = [ "@ERR445333.1.2 1 length=30",
        "AACCCCTAACCCCTAACCCTAACCCTAACC",
        "+ERR445333.1.2 1 length=30",
        "EFFFGF?EFGA???DC>:?A##########" ]
  resultLines = convertFastqItem2FastaItem(oneFastqItem, 10)
  assert(2, len(resultLines))
  assert('>ERR445333.1.2 1 length=20' == resultLines[0])
  assert(20, len(resultLines[1])) # last 10 chars cut 

def test_trim_all_bad():
  oneFastqItem = [ "@ERR445333.1.2 1 length=30",
        "AACCCCTAACCCCTAACCCTAACCCTAACC",
        "+ERR445333.1.2 1 length=30",
        "##############################" ]
  resultLines = convertFastqItem2FastaItem(oneFastqItem, 10)
  assert("TOO LOW SEQUENCE QUALITY", resultLines)

def test_trim1():
  oneFastqItem = [ "@ERR445333.1.2 1 length=30",
        "AACCCCTAACCCCTAACCCTAACCCTAACC",
        "+ERR445333.1.2 1 length=30",
        "EFFFGF?EFGA???DC>:?AEA########" ]
  resultLines = convertFastqItem2FastaItem(oneFastqItem, 20, True)
  assert(2, len(resultLines))
  assert('>ERR445333.1.2 1 length=22' == resultLines[0])
  assert(22, len(resultLines[1])) # nothing should be cut

def test_notrim1():
  oneFastqItem = [ "@ERR445333.1.2 1 length=30",
        "AACCCCTAACCCCTAACCCTAACCCTAACC",
        "+ERR445333.1.2 1 length=30",
        "EFFFGF?EFGA???DC>:?AEA########" ]
  resultLines = convertFastqItem2FastaItem(oneFastqItem, 20, False)
  assert(2, len(resultLines))
  assert('>ERR445333.1.2 1 length=30' == resultLines[0])
  assert(30, len(resultLines[1])) # nothing should be cut
  
