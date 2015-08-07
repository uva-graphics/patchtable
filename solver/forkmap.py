
"""
forkmap -- Forking map(), uses all processors by default.

Connelly Barnes 2008, public domain.  Based on forkmap by Kirk Strauser, rewritten and optimized.  Version 1.0.2.
"""

import os, mmap, marshal, struct, cPickle
import ctypes, ctypes.util
import time, traceback

builtin_map = map

def nprocessors():
  try:
    try:
      try:
        # Multiprocessing
        import multiprocessing
        return multiprocessing.cpu_count()
      except:
        # Mac OS
        libc=ctypes.cdll.LoadLibrary(ctypes.util.find_library('libc'))
        v=ctypes.c_int(0)
        size=ctypes.c_size_t(ctypes.sizeof(v))
        libc.sysctlbyname('hw.ncpu', ctypes.c_voidp(ctypes.addressof(v)), ctypes.addressof(size), None, 0)
        return v.value
    except:
      # Cygwin (Windows) and Linuxes
      # Could try sysconf(_SC_NPROCESSORS_ONLN) (LSB) next.  Instead, count processors in cpuinfo.
      s = open('/proc/cpuinfo', 'r').read()
      return s.replace(' ', '').replace('\t', '').count('processor:')
  except:
    return 1

nproc = nprocessors()

def map(f, *a, **kw):
  """
  forkmap.map(..., n=nprocessors), same as map(...).

  n must be a keyword arg; default n is number of physical processors.
  """
  def writeobj(pipe, obj):
    try:
      s = marshal.dumps(obj)
      s = struct.pack('i', len(s)) + s
    except:
      s = cPickle.dumps(obj)
      s = struct.pack('i', -len(s)) + s
    os.write(pipe, s)

  def readobj(pipe):
    n = struct.unpack('i', os.read(pipe, 4))[0]
    s = ''
    an = abs(n)
    while len(s) < an:
      s += os.read(pipe, min(65536, an-len(s)))
    if n > 0:
      return marshal.loads(s)
    else:
      return cPickle.loads(s)

  n = kw.get('n', nproc)
  if n == 1:
    return builtin_map(f, *a)

  if len(a) == 1:
    L = a[0]
  else:
    L = zip(*a)
  try:
    len(L)
  except TypeError:
    L = list(L)
  n = min(n, len(L))

  ans = [None] * len(L)
  pipes = [os.pipe() for i in range(n-1)]

  for i in range(n):
    if i < n-1 and not os.fork():
      # Child, and not last processor
      try:
        try:
          if len(a) == 1:
            obj = builtin_map(f, L[i*len(L)//n:(i+1)*len(L)//n])
          else:
            obj = [f(*x) for x in L[i*len(L)//n:(i+1)*len(L)//n]]
        except Exception, obj:
          pass
        writeobj(pipes[i][1], obj)
      except:
        traceback.print_exc()
      finally:
        os._exit(0)
    elif i == n-1:
      # Parent fork, and last processor
      try:
        if len(a) == 1:
          ans[i*len(L)//n:] = builtin_map(f, L[i*len(L)//n:])
        else:
          ans[i*len(L)//n:] = [f(*x) for x in L[i*len(L)//n:]]
        for k in range(n-1):
          obj = readobj(pipes[k][0])
          if isinstance(obj, Exception):
            raise obj
          ans[k*len(L)//n:(k+1)*len(L)//n] = obj
      finally:
        for j in range(n-1):
          os.close(pipes[j][0])
          os.close(pipes[j][1])
          os.wait()
  return ans

def bench():
  print 'Benchmark:\n'
  def timefunc(F):
    start = time.time()
    F()
    return time.time() - start
  def f1():
    return builtin_map(lambda x: pow(x,10**1000,10**9), range(10**3))
  def g1():
    return map(lambda x: pow(x,10**1000,10**9), range(10**3))
  def f2():
    return builtin_map(lambda x: x**2, range(10**6))
  def g2():
    return map(lambda x: x**2, range(10**6))
  import timeit
  print 'Expensive operation, 10**3 items:'
  print 'map         (1 processor): ', timefunc(f1), 's'
  print 'forkmap.map (%d processors):' % nproc, timefunc(g1), 's'
  print
  print 'Cheap operation, 10**6 items:'
  print 'map         (1 processor): ', timefunc(f2), 's'
  print 'forkmap.map (%d processors):' % nproc, timefunc(g2), 's'

def test():
  print 'Testing:'
  assert [x**2 for x in range(10**4)] == map(lambda x: x**2, range(10**4))
  assert [x**2 for x in range(10**4)] == map(lambda x: x**2, range(10**4), n=10)
  assert [x**2 for x in range(10**4)] == map(lambda x: x**2, range(10**4), n=1)
  assert [(x**2,) for x in range(10**3,10**4)] == map(lambda x: (x**2,), range(10**3,10**4))
  assert [(x**2,) for x in range(10**3,10**4)] == map(lambda x: (x**2,), range(10**3,10**4), n=10)
  assert [(x**2,) for x in range(10**3,10**4)] == map(lambda x: (x**2,), range(10**3,10**4), n=1)
  assert builtin_map(lambda x,y:x+2*y, range(100),range(0,200,2)) == map(lambda x,y:x+2*y, range(100),range(0,200,2))
  assert builtin_map(lambda x,y:x+2*y, range(100),range(0,200,2)) == map(lambda x,y:x+2*y, range(100),range(0,200,2), n=10)
  assert builtin_map(lambda x,y:x+2*y, range(100),range(0,200,2)) == map(lambda x,y:x+2*y, range(100),range(0,200,2), n=2)
  # Some Windows (Cygwin) boxes can't fork more than about 15 times, so only test to n=15
  for n in range(1, 15):
    assert [x**3 for x in range(200)] == map(lambda x: x**3, range(200), n=n)
  def f(n):
    if n == 1:
      raise KeyError
  def check_raises(func, exc):
    e = None
    try:
      func()
    except Exception, e:
      pass
    if not isinstance(e, exc):
      raise ValueError('function did not raise specified error')

  check_raises(lambda: map(f, [1, 0], n=2), KeyError)
  check_raises(lambda: map(f, [0, 1], n=2), KeyError)
  check_raises(lambda: map(f, [1, 0, 0], n=3), KeyError)
  check_raises(lambda: map(f, [0, 1, 0], n=3), KeyError)
  check_raises(lambda: map(f, [0, 0, 1], n=3), KeyError)
  print 'forkmap.map: OK'

if __name__ == '__main__':
  test()
  bench()
