"""
Copyright (c) 2005 Pennsylvania State University

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Some icons found in Galaxy are from the Silk Icons set, available under
the Creative Commons Attribution 2.5 License, from:

http://www.famfamfam.com/lab/icons/silk/

"""


#!/usr/bin/env python
#Guruprasad Ananda
#modifications to support python3 by sebastien weyn
"""
This tool provides the UNIX "join" functionality.
"""
import sys, os, tempfile, subprocess

def stop_err(msg):
	sys.stderr.write(msg)
	sys.exit()

def main():
	infile1 = sys.argv[1]
	infile2 = sys.argv[2]
	field1 = int(sys.argv[3])
	field2 = int(sys.argv[4])
	mode =sys.argv[5]
	outfile = sys.argv[6]
	
	tmpfile1 = tempfile.NamedTemporaryFile()
	tmpfile2 = tempfile.NamedTemporaryFile()
	
	try:
		#Sort the two files based on specified fields
		os.system("sort -t '	' -k %d,%d -o %s %s" %(field1, field1, tmpfile1.name, infile1))
		os.system("sort -t '	' -k %d,%d -o %s %s" %(field2, field2, tmpfile2.name, infile2))
	except:
		stop_err( 'Initialization error -> %s' %str(exc) )
		
	option = ""
	with open(tmpfile1.name, 'r') as infile:
		for line in infile:
			line = line.strip()
			if line:
				elems = line.split('\t')
				for j in range(1,len(elems)+1):
					if j == 1:
						option = "1.1"
					else:
						option = option + ",1." + str(j) 
			break
	
	#check if join has --version option. BSD join doens't have this option, while GNU join does. 
	#The return value in the latter case will be 0, and non-zero in the latter case.
	ret = subprocess.run('join --version 2>/dev/null', shell=True) 
	# check if we are a version later than 7 of join. If so, we want to skip
	# checking the order since join will raise an error with duplicated items in
	# the two files being joined.
	if ret.returncode == 0: 
		cl = subprocess.Popen(["join", "--version"], stdout=subprocess.PIPE)
		(stdout, _) = cl.communicate()
		version_line = str(stdout).split("\\n")[0]
		(version, _) = version_line.split()[-1].split(".")
		if int(version) >= 7:
			flags = "--nocheck-order"
		else:
			flags = ""
	else:
		flags = ""

	if mode == "V":
		cmdline = "join %s -t '	' -v 1 -o %s -1 %d -2 %d %s %s > %s" %(flags, option, field1, field2, tmpfile1.name, tmpfile2.name, outfile)
	else:
		cmdline = "join %s -t '	' -o %s -1 %d -2 %d %s %s > %s" %(flags, option, field1, field2, tmpfile1.name, tmpfile2.name, outfile)
	
	try:
		os.system(cmdline) 
	except:
		stop_err('Error joining the two datasets -> %s' %str(exj))
	   
if __name__ == "__main__":
	main()
