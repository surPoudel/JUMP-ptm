import re


def read_header_pep_xml(lines_list,pep_xml_dictionary):
    for index, line in enumerate(lines_list):
        if '</search_summary>' in line.strip():
            pep_xml_dictionary["header"].append(line)
            break
        else:
            pep_xml_dictionary["header"].append(line)
    return index


def find_scan_pep_xml(line, fraction):
    pattern = '<spectrum_query spectrum="('+fraction+'\.\d+\.\d+\.\d+)".+precursor_neutral_mass="(\d+\.\d+)"'
    allPat=re.match(pattern, line.strip())
    x = allPat.group(1)
    return x #x is the spectrum



def read_scan_store(lines_list):
    pep_xml_dictionary = {"header":[]}
    #store headers
    index = read_header_pep_xml(lines_list,pep_xml_dictionary)
    #store scans

    #</spectrum_query>\n
    spectrum = None
    for x in lines_list[index:]:

        if '<spectrum_query spectrum' in x.strip(): #scan found
            spectrum = find_scan_pep_xml(x, "w001")
            pep_xml_dictionary[spectrum] = []
        if ('</search_summary>' not in x.strip()) and (spectrum!=None):
    #         print (spectrum,"\t",x.strip())
            pep_xml_dictionary[spectrum].append(x)
    
    #add tail last 2 lines ' </msms_run_summary>\n','</msms_pipeline_analysis>\n' to last spectrum key
    
    pep_xml_dictionary[spectrum]+lines_list[-2:]
    
    return pep_xml_dictionary


