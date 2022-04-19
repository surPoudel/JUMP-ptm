from itertools import product


def filler(word, from_char, to_char):
    options = [(c,) if c != from_char else (from_char, to_char) for c in word]
    return (''.join(o) for o in product(*options))



def merge_two_dicts(x, y):
    z = x.copy()   # start with keys and values of x
    z.update(y)    # modifies z with keys and values of y
    return z




#a function to parse tagFile. This will create a dictionary that has the spectrum as the key. 
#Each spectrum dictionary has another dictionary as the values that has Tag as key and its score as values

#w010.1000.1.2.tag
#1279.68799  ENGDAST 904.419366424669
def scan_tag_dictionary(tag_file): #change key to outfile style with precursor info
    tag_score_list = []
    scan_tag = {}
    jump_comet_link = {} #links comet with jump
    tag_rank = 0
    with open(tag_file,"r") as f:
        cnt = 0
        for line in f:
            tag_rank+=1
            # print (cnt)
            if '.tag' in line:
                tag_rank = 0
                temp_line = line.split(".")
                scan = int(temp_line[1])
                scan5digit = str('%05d' % int(scan)) 
                #this is comet key
                comet_key = temp_line[0]+"."+scan5digit+"."+scan5digit+"."+temp_line[-2]
                jump_key = ".".join(temp_line[0:-1])
                tag_score_dict = {} #tag score dictionary
                link_key = ""
                
            else:
                
                # print (tag_rank)
                tag_info = line.strip().split("\t")
                neutral_mass = tag_info[0] #MH -- hydrogen is added to neutral mass
                link_key = comet_key+"_"+neutral_mass #we are adding w010.00500.00500.3 with neutral mass (MH)
                tag = tag_info[1]
                rtFlank = tag_info[2]
                Escore = tag_info[3] #use second score in .tags file
#                 if key == "w010.06508.06508.3":
#                     print (tag_info)
#                     print (rtFlank,"\t",Escore)
                
                
                if "@" in tag:
                    tag1 = tag.replace('M@','M')
                    tag=tag1
                if "I" in tag or "L" in tag:
                    leu_pep = tag.replace('I','L')
                    pep_list = list(filler(leu_pep,'L','I'))

                    for tags in pep_list:
                        tag_score_dict[tags] = [Escore,rtFlank,tag_rank]
#                         tag_score_dict[tags[::-1]] = [Escore,rtFlank] #jump already reverses the tag
                        tag_score_list.append(float(Escore))
                else:
                    tag_score_dict[tag] = [Escore,rtFlank,tag_rank]
#                     tag_score_dict[tag[::-1]] = [Escore,rtFlank]
                    tag_score_list.append(float(Escore))
            scan_tag[jump_key] = tag_score_dict
            jump_comet_link[link_key] = jump_key

            # if jump_key not in scan_tag:
            #     scan_tag[jump_key] = tag_score_dict
            # else:
            #     #update tag_score_dict
            #     scan_tag[jump_key] = merge_two_dicts(scan_tag[jump_key], tag_score_dict)

    return scan_tag, tag_score_list, jump_comet_link

#dict1[key] = [val]
#all_matched_tags, spectrum, [key,ionType,position]
#key,ionType,rank,tag_index, start_position, end_position,tag_score,rtflank_tag_rank
def report_tag(spectrum_tag_dict, out_file):
    with open(out_file,"w") as f:
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("spectrum","tag","ionType","rank","tag_index","tag_start_position","tag_end_position","total_tags","tag_score","tag_rank"))
        for key in spectrum_tag_dict.keys():
            alltags = spectrum_tag_dict[key]
            for val in alltags:
                # val_split = val[0].split("\t")
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key,val[0],val[1],val[2],val[3],val[4],val[5],len(alltags),val[6],val[7]))

