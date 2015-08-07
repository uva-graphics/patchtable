import time
import subprocess
model = subprocess.check_output(['bash', '-c', 'cat /proc/cpuinfo | grep "model name" | head -n 1'])
host = subprocess.check_output(['bash', '-c', 'hostname'])
s = model.strip() + '|' + host.strip()

with open('loghost.txt', 'at') as f:
    print >> f, s
time.sleep(30.0)
