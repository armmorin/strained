from pathlib import Path
import re

cwd = Path.cwd()
list_of_files = list(cwd.glob('perqueue*err'))

def search_file(file):
    for file in list_of_files: 
        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if re.search(r'CalculationFailed', line):
                    print(file, line)
                    return

search_file(list_of_files)
        