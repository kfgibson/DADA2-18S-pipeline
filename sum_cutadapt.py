cutadapt_sum = "./cutadapt_stats.txt"

print("Sample name\treads passing\tbases passing")

with open(cutadapt_sum, "r") as in_file:
	for line in in_file:
		line = line.rstrip()
		if "Command" in line:
			sample = line.split()[-1]
		if "filters" in line:
			reads_passing = line.split("(")[-1][:-1]
		if "filtered" in line:
			bases_passing = line.split("(")[-1][:-1]
			print(sample + "\t" + reads_passing + "\t" + bases_passing)
