import stat
import os

# if on the host, both will print False
# but if in the container, the first will print True
mode = os.fstat(0).st_mode
print(stat.S_ISFIFO(mode))
print(stat.S_ISREG(mode))
