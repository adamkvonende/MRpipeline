<br>

### Harmonising on chr:pos

<br>
Currently, the pipeline will only harmonise on rsid. However, there
is a workaround if you wish to harmonise on chr:pos. You can do this by
creating an 'rsid' column in the exposure and outcome data sets that
includes chr:pos (e.g. chr11:14335353). The reason this works is because
the harmonisation function simply matches strings from one data set to
another, so does not care if it is 'rs14343' or 'chr11:14335353' or
'flyingmonkey1234'. There is one catch: if you want the function to
perform clumping or create an LD matrix for you, these will not work
using chr:pos because these functions rely on rsids. However, there is
now a workaround for this as well: we have added an experimental
argument ('map_rsids=T') that will map chr:pos from harmonised SNPs to
the rsid, so that the function can then perform downstream clumping,
etc,. You will need two additional packages for this:
SNPlocs.Hsapiens.dbSNP144.GRCh37 and SNPlocs.Hsapiens.dbSNP155.GRCh38.
If you are harmonizing on chr:pos, make sure your exposure and outcome
data sets are on the same build!