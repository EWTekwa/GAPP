infile = open("results.tab")
outfile = open("out.tab", "w")
for line in infile:
    col = i + "\t" + line
    columns = col.split("\t")
    print(columns)
    print(col)
    outfile.write("\t".join(columns)+"\n")