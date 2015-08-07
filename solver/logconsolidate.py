import sys

s = open('loghost.txt','rt').read()
L = s.strip().split('\n')
host = []
model = []
for line in L:
  suffix = '.local'
  h = line.split('|')[1]
  if h.endswith(suffix):
    h = h[:len(h)-len(suffix)]
  host.append(h)
  m = line.split('|')[0]
  m = m[m.index('CPU')+3:].strip().split()[0]
  model.append(m)

d = {}
dset = {}
for i in range(len(host)):
  d.setdefault(model[i], [])
  d[model[i]].append(host[i])
  dset.setdefault(model[i], set())
  dset[model[i]].add(host[i])

for key in dset:
  dset[key] = sorted(dset[key])

for key in sorted(d.keys()):
  print >> sys.stderr, '%-10s'%key, 'unique nodes:', '%-02d'%len(set(d[key])), '  count in log file:', len(d[key])

print dset

