import os
import subprocess

import sys

src_dir = os.getenv('SRC_DIR')

test_command = "runUnitTests"
print("Test command = %s" % test_command)
command_result = subprocess.call(test_command, shell=True)
print("ran test command")
sys.exit(command_result)