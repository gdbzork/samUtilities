from os.path import splitext
import logging
import pysam

class BetterHits:

  LOGNAME = "betterHits"
  LOGLEV = logging.DEBUG

  def __init__(self,log=None):
    self.log = log if log is not None else self.configureLogging()
    self.results = None

  def configureLogging(self):
    log = logging.getLogger(self.LOGNAME)
    log.setLevel(self.LOGLEV)
    fmt = logging.Formatter(fmt="%(asctime)s %(levelname)s %(message)s",
                            datefmt="[%Y-%m-%d %H:%M:%S]")
    hdlr = logging.StreamHandler()
    hdlr.setFormatter(fmt)
    log.addHandler(hdlr)
    return log

  @staticmethod
  def openFile(fn):
    # test suffix: is it BAM or SAM?
    suff = splitext(fn)[1]
    if suff == ".bam":
      fd = pysam.AlignmentFile(fn,"rb")
    elif suff == ".sam":
      fd = pysam.AlignmentFile(fn,"r")
    else:
      raise OSError("unknown suffix '%s' for '%s'" % (suff,fn))
    return fd

  @staticmethod
  def loadGTF(fn,regions):
    with open(fn,'r') as fd:
      for line in fd:
        flds = line.split("\t")
        regions.append((flds[0],int(flds[3]),int(flds[4])))
    return regions

  def loadHits(self,fn):
    fd = self.openFile(fn)
    hits = {}
    count = 0
    primary = 0
    for read in fd.fetch(until_eof=True):
      count += 1
      if read.is_secondary:
        continue
      primary += 1
      score = read.get_tag("AS")
      hits[read.query_name] = score
    fd.close()
    self.log.info("loaded %d reference sequences (%d primary)" % (count,primary))
    return hits

  @staticmethod
  def overlaps(read,intervals):
    olap = False
    for (chrom,left,right) in intervals:
      if chrom == read.reference_name:
        olapBP = read.get_overlap(left,right)
        if olapBP > 0:
          olap = True
          break
    return olap

  def getResults(self):
    return self.results

  def dumpResults(self,dest):
    for (name,ascore,rscore) in self.results:
      dest.write("%s\t%d\t%d\n" % (name,ascore,rscore))

  def compare(self,referenceFN,alternativeFN,gtfFN):
    ref = self.loadHits(referenceFN)
    altFD = self.openFile(alternativeFN)
    intervals = []
    if gtfFN is not None:
      self.loadGTF(gtfFN,intervals)
    acount = 0
    aproblem = 0
    foundInAlt = set()
    results = []
    for alt in altFD.fetch():
      acount += 1
      aname = alt.query_name
      if aname not in ref:
        continue
      if self.overlaps(alt,intervals):
        continue
      foundInAlt.add(aname)
      ascore = alt.get_tag("AS")
      if ascore > ref[aname]:
        aproblem += 1
        results.append((aname,ascore,ref[aname]))
    self.results = results
    self.log.info("checked %d alt sequences (%d problems, %d foundInAlt)\n" % (acount,aproblem,len(foundInAlt)))
    return len(results)
