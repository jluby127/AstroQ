import numpy as np
import pandas as pd
import sys
import argparse

font = '\\fontsize{\\fsize}{14}\selectfont\parbox[c]{\linewidth}{'
bold = '\\textbf{'
templateFile = 'CoversheetTemplate.txt'
path2requests = '/Users/jack/Desktop/'

def writeRequestReport(semester, program):
    data = pd.read_csv(path2requests + str(semester) + "_" + str(program) + "_Requests.csv")
    totaltime = np.round(np.sum(data['Total Time for Target [hours].1']),2)

    tableLines = []
    for x in range(len(data)):
        if data['Request Feasibility'][x]:
            feasible = "Yes"
        else:
            feasible = 'No'
        cell1 = font + "\centering " + str(data['Target Name'][x]) + "} & "
        cell2 = font + "\centering " + str(data['Nominal Exposure Time [s]'][x]) + "} & "
        cell3a = '\hfill \\\ \hfill \\\ \centering On \\textbf{' + str(data['# of Nights Per Semester'][x]) + '} unique nights \\newline '
        cell3b = 'at \\textbf{' + str(data['Minimum Inter-Night Cadence'][x]) + '} day minimum inter-night cadence \\newline '
        cell3c = 'observe star \\textbf{' + str(data['# Visits per Night'][x]) + '} visit(s) \\newline '
        cell3d = 'at \\textbf{' + str(data['Minimum Intra-Night Cadence'][x]) + '} hr minimum intra-night cadence \\newline '
        cell3e = 'and take \\textbf{' + str(data['# of Exposures per Visit'][x]) + '} exposure(s) each time \\\ \hfill \\\} &  '
        cell3 = font + cell3a + cell3b + cell3c + cell3d + cell3e
        cell4 = font + "\centering " + str(data['Total Time for Target [hours].1'][x]) + "} & "
        cell5 = font + "\centering " + str(feasible) + "} \\\ " # $^a$
        totalline = cell1 + cell2 + cell3 + cell4 + cell5
        tableLines.append(totalline)

    file_original = open(templateFile, 'r')
    Lines = file_original.readlines()
    file_original.close()

    file_new = open(path2requests + str(semester) + "_" +str(program) + "RequestReport.txt", 'w')
    for line in Lines:
        if line.strip() == 'Semester ID:':
            file_new.write('Semester ID: ' + bold + str(semester) + ' ' + str(program) + '} \n')
        elif line.strip() == 'Total Time Required to Complete Program (hr):':
            file_new.write('Total Time Required to Complete Program (hr): ' + bold + str(totaltime) + '} \n')
        else:
            file_new.write(line)

    for i in range(len(tableLines)):
        file_new.write(tableLines[i] + "\n")
        file_new.write("\hline \n")

    file_new.write('\hline \n')
    file_new.write('\end{longtable} \n')
    file_new.write('\end{document} \n')
    file_new.close()

parser = argparse.ArgumentParser(description='Generate request report tex file for TAC')
parser.add_argument('-s','--semester', help='The 4 digit year and 1 character semester letter')
parser.add_argument('-p','--program', help='The 4 character program code')
args = parser.parse_args()

print("Generating Request Report for " + str(str(args.semester) + "_" + str(args.program)) + ".")
writeRequestReport(args.semester, args.program)
print("Done.")
