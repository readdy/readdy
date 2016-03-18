import os
import subprocess

import sys

src_dir = os.getenv('SRC_DIR')

test_command = "cd %s/build; make runUnitTests; make test ARGS=\"-V\"; cd -" % src_dir
print("Test command = %s" % test_command)
command_result = subprocess.call(test_command, shell=True)
print("ran test command")
sys.exit(command_result)