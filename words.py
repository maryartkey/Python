textfile = open("sonnets.txt", 'r') 
textoutput = open("title.txt", 'w')

line = textfile.readline()
while line:
	if line == '\n':
		line = textfile.readline()
		textoutput.write(line)
		print(line)
		line = textfile.readline()
		line = textfile.readline()
	line = textfile.readline()
textfile.close()
