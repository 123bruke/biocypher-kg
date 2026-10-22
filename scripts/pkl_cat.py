import gzip
import pickle
import sys

# Check if a filename was passed as an argument
if len(sys.argv) < 2:
    print("Usage: python3 pkl_cat.py <file_name.pkl>")
    sys.exit(1)

file_path = sys.argv[1]

GZIP_MAGIC = b"\x1f\x8b"

try:
    with open(file_path, 'rb') as f:
        is_gzipped = f.read(2) == GZIP_MAGIC

    opener = gzip.open if is_gzipped else open
    with opener(file_path, 'rb') as f:
        data = pickle.load(f)
        print("File:", file_path)
        print("Object type:", type(data))
        print("-" * 20)
        print("Content:")
        for tissue in data:
            print((tissue, data.get(tissue)))
except FileNotFoundError:
    print(f"Error: File not found at '{file_path}'")
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
