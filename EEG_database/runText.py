import subprocess
import shlex

def run_script_from_file(filename):
    # Hardcoded script path (absolute or relative)
    script_path = r"C:\Users\samckenzie\PycharmProjects\deid_edf\main.py"

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines, start=1):
                line = line.strip()
                if not line:
                    continue

                # Split line by whitespace into parts
                parts = shlex.split(line)


                # We expect at least 3 parts: script_name, arg1, arg2
                if len(parts) < 3:
                    print(f"Line {i} is malformed (less than 3 parts): {line}")
                    continue

                # Extract the two arguments only (ignore first column)
                arg1 = parts[1]
                arg2 = parts[2]

                # Construct command: python script_path arg1 arg2
                cmd = ["python", script_path, arg1, arg2]

                print(f"Running line {i}: {cmd}")
                try:
                    subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error on line {i}: {e}")
    except FileNotFoundError:
        print(f"The file '{filename}' was not found.")


if __name__ == "__main__":
    badfils = r"N:\Research-Studies\Study 21-044 HRDI BM\Code\fix_bad_fils7.txt"
    run_script_from_file(badfils)