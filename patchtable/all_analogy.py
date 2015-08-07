
import os

cmd = './patchtable analogy image_analogies/%s/a.png image_analogies/%s/ap.png b.png out/%s.png -threads 8 -coherence_spatial 80 -upweight_coherence 1'

L = """
engraving_brick
engraving_moon
engraving_nast
engraving_nast_large
engraving_squire
engraving_waterfall
freud
manet
pastel
stipple
stipple_chair
watercolor
""".strip().split()

L = """
penandink
watercolor2
""".strip().split()

def system(s):
  print s
  os.system(s)

for dir in L: #['stipple_chair']: #L:
  system(cmd % (dir, dir, dir))

