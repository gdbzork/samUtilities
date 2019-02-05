samUtilities
============

Miscellaneous bits and pieces for messing with SAM- and BAM-format files.

Python scripts
--------------

:isPrimaryAlwaysBest: test whether the alignment marked as "primary" is
                      actually among the highest-scoring hits
:keepOffHits: keep hits that are imperfect (NM > 0 or cigar doesn't
              match '^\d+M$') for later analysis

:reportBetterHits: report hits in second file that score higher than
                   first.
:samflags: decode SAM flags field
