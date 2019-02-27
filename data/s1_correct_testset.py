import sys




testfile=sys.argv[1]
test_r1_file=sys.argv[2]
outfile=sys.argv[3]

lines_split= open(testfile).readlines()
lines1= open(test_r1_file).readlines()


for line in lines1[1:]:
    key=int(line.split(",")[0])+1
    lines_split[key]=line


with open(outfile, "w") as f:
    f.write("".join(lines_split))



