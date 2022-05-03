fs = 2000
n_mins = 20

csv_in_path = sys.argv[1] # csv path passed as argument to the script
head,csv_fname = os.path.split(csv_path)
basename,_ = os.path.splitext(csv_fname)
csv_out_path = os.path.join(head,f"{basename}_{n_mins}.csv")

with open(csv_in_path,'r') as f, open(csv_out_path,'w') as fnew:
    for i in range(fs*n_mins):
        fnew.write(f.readline())


