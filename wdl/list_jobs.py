#!/usr/bin/env python3
import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

os.system("cromshell list -c -u 2>&1 | tail -n+8")
