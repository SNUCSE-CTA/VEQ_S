import sys

import libVEQ

if len(sys.argv) != 5:
    print("Invalid arguments")

libVEQ.run(
    sys.argv[1:],
)
