import os
import subprocess
import sys

print("***"*15 + " running c++ tests " + "***"*15)

test_command = "runUnitTests"
print("Test command = %s" % test_command)
command_result = subprocess.call(test_command, shell=True)
print("ran test command")
if command_result != 0:
    sys.exit(command_result)

print("***"*15 + " running singlecpu c++ tests " + "***"*15)

test_command = "runUnitTests_singlecpu"
print("Test command = %s" % test_command)
command_result = subprocess.call(test_command, shell=True)
print("ran test command")
if command_result != 0:
    sys.exit(command_result)


print("***"*15 + " running python tests " + "***"*15)

src_dir = os.getenv('SRC_DIR')
test_pkg = 'readdy'
nose_run = "nosetests {test_pkg} -vv" \
           " --with-doctest --doctest-options=+NORMALIZE_WHITESPACE,+ELLIPSIS" \
    .format(test_pkg=test_pkg).split(' ')
command_result = subprocess.call(nose_run)
sys.exit(command_result)
