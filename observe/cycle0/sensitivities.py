
from math import *

# NB for best part of band

###################################

# Set these variables
# Start python or casa
# execfile('sensitivities.py')
# returns data needed to update Cycle0.html

# 12 hr on-source (#18 hr total track)
# uJy/beam

# 320 MHz total b/w, ~300 useful for continuum
Lbw = 320           # total bandwidth
Lnbands = 10        # Number of sub-bands in total
LwithLo = 6.        # uJy/bm with Lo
LnoLo   = 12.       # uJy/bm with Lo

# 512 MHz b/w, 500 useful for continuum
Cbw = 512         # total bandwidth
Cnbands = 4.      # Number of sub-bands in total
CwithLo = 7.      # uJy/bm with Lo
CnoLo   = 13.     # uJy/bm with Lo

########################################
f = sqrt(2.)

# Full bandwidth

print '320 MHz total bandwidth'
print 'L-band 12-hr, with Lo, uJy/b', '%4.0f' % (LwithLo)
print 'L-band 12-hr,   no Lo, uJy/b', '%4.0f' % (LnoLo)
#print 'L-band 6-hr,  with Lo, uJy/b', '%4.0f' % (LwithLo*f)
#print 'L-band 6-hr,    no Lo, uJy/b', '%4.0f' % (LnoLo*f)

print '512 MHz total bandwidth'
print 'C-band 12-hr, with Lo, uJy/b', '%4.0f' % (CwithLo)
print 'C-band 12-hr,   no Lo, uJy/b', '%4.0f' % (CnoLo)
#print 'C-band 6-hr,  with Lo, uJy/b', '%4.0f' % (CwithLo*f)
#print 'C-band 6-hr,    no Lo, uJy/b', '%4.0f' % (CnoLo*f)

# Spectral, L-band with Lo
b0 = Lbw/Lnbands

print '\nL-band, with Lo, 12 hr'
print 'spw MHz mJy/bm/ch 1000K'

for i in range(int(3+log(b0)/log(2))):
    spw = b0/2.**i
    s = LwithLo*sqrt(320./spw)*sqrt(512.)/1000.
    print '%6.2f   %6.2f   %5.1f' % (spw,  s,   s*15.)

# Spectral, L-band no Lo
b0 = Lbw/Lnbands

print '\nL-band, with Lo, 12 hr'
print 'spw MHz mJy/bm/ch 1000K'

for i in range(int(3+log(b0)/log(2))):    
    spw = b0/2.**i
    s = LnoLo*sqrt(320./spw)*sqrt(512.)/1000.
    print '%6.2f   %6.2f   %5.1f' % (spw,  s,   s*15.)

# Spectral, C-band with Lo
b0 = Cbw/Cnbands

print '\nC-band, with Lo, 12 hr'
print 'spw MHz mJy/bm/ch 1000K'

for i in range(int(3+log(b0)/log(2))):
    spw = b0/2.**i
    s = CwithLo*sqrt(320./spw)*sqrt(512.)/1000.
    print '%6.2f   %6.2f   %5.1f' % (spw,  s,   s*15.)

# Spectral, C-band no Lo
b0 = Cbw/Cnbands

print '\nC-band, with Lo, 12 hr'
print 'spw MHz mJy/bm/ch 1000K'

for i in range(int(3+log(b0)/log(2))):
    spw = b0/2.**i
    s = CnoLo*sqrt(320./spw)*sqrt(512.)/1000.
    print '%6.2f   %6.2f   %5.1f' % (spw,  s,   s*15.)
