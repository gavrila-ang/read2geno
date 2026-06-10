'''Settings for formatting print statements'''

# ANSI color codes
R = '\033[91m' # red
B = '\033[94m' # blue
G = '\033[92m' # green
D = '\033[0m'  # reset to default color (D = 'Done')


# static symbols - can't be modified
PASS = '✅'  # check mark
FAIL = '❌' # cross mark
WARN = '⚠️' # warning symbol
FILE = '📝' # memo (for writing files)

# generic check & cross marks
V = '✔'
X = '✘'

# set colors for check/cross marks
V = f'{B}{V}{D}'
X = f'{R}{X}{D}'