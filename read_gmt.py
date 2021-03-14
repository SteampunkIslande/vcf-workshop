import sys

for line in sys.stdin:
    sys.stdout.write(f"{len(line.split(','))} {line.split(',')[0]}\n")

