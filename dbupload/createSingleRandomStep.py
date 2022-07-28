import AccessDb
import sys

a=AccessDb.AccessDb("TI_27-JAN-2010_2")
a.Download()
a.DelayDetIds("<name of text files with delays>","<name of XML file for the partition to be produced with PLL and FED coarse and fine>")
a.setUploadFeds(True)
a.setUploadFecs(True)
a.Upload()

b=AccessDb.AccessDb("TO_30-JUN-2009_1")
b.Download()
a.DelayDetIds("<name of text files with delays>","<name of XML file for the partition to be produced with PLL and FED coarse and fine>")
b.setUploadFeds(True)
b.setUploadFecs(True)
b.Upload()

c=AccessDb.AccessDb("TP_09-JUN-2009_1")
c.Download()
a.DelayDetIds("<name of text files with delays>","<name of XML file for the partition to be produced with PLL and FED coarse and fine>")
c.setUploadFeds(True)
c.setUploadFecs(True)
c.Upload()

d=AccessDb.AccessDb("TM_09-JUN-2009_1")
d.Download()
a.DelayDetIds("<name of text files with delays>","<name of XML file for the partition to be produced with PLL and FED coarse and fine>")
d.setUploadFeds(True)
d.setUploadFecs(True)
d.Upload()
