
def test_read_fastq():
  in_data = open('../sample_data/sample1.fastq', 'r')
  counter = 0
  for line in in_data:
    counter += 1
  assert(200, counter)

