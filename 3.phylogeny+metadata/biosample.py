#!/usr/bin/python

import re, csv

def unique_attributes(file_path):
    unique_strings = set()

    with open(file_path, 'r') as file:
        content = file.read()
        matches = re.findall(r'/.*=', content)
        unique_strings.update(matches)

    attributes = ['Accession', 'Number'] + list(unique_strings)
    attributes = [attr.strip('/=') for attr in attributes]

    return attributes

def manipulate_table(inputname, outputname):
    #reading file
    input = open(inputname, 'r')
    sample_dict = {}

    with open(outputname, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=unique_attributes(inputname), restval='NA', extrasaction='raise')
        writer.writeheader()
        #read blocks
        for line in input:
            if process_line(sample_dict, line):
                #dictionary in tabelle schreiben
                writer.writerow(sample_dict)
                #dict leeren
                sample_dict.clear()

    input.close()

def process_line(sample_dict, line):
    if line[0].isdigit(): #if starts with number: extract number
        number = line.split(':')[0]
        sample_dict["Number"] = number
    elif line.strip().startswith(r'/'): #if starts with /attribute: extract attribute
        line = line.strip().strip('/')
        attribute_pair = line.split('=')
        sample_dict[attribute_pair[0]] = attribute_pair[1].strip('"')
        if sample_dict.get(attribute_pair[0]) == "missing":
            sample_dict[attribute_pair[0]] = 'NA'
    elif line.startswith('Accession'): #if starts with Accession: extract accession
        accession = line.split()[1]
        sample_dict['Accession'] = accession
    elif line in ['\n', '\r\n']: #if empty line: block ends
        return True
    else:
        pass #line has no interesting information
    return False
    

def main():
    #print(unique_attributes('biosample_results.txt')[:5])
    manipulate_table('biosample_results.txt', 'biosample_table.csv')

if __name__ == "__main__":
    main()