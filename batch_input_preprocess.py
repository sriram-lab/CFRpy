def batch_input_preprocess(data_dir, prefix_name, medium):
    '''Extracts prefixes and data series from single-cell datasets
    data_dir (str): file path to directory containing files to process
    prefix_name (str): prefix to be applied to file names
    medium (str): name of extracellular metabolite source ('DMEMF12' or 'KSOM')'''

    data_dir_files = sorted([f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f))]) # get files in data directory

    data_series = []
    prefix_series = []

    for file_name in data_dir_files:

        # create data_series
        data_name = file_name.split('_') # convert file_name into list of substrings sep by _

        if len(data_name) > 1 and 'genes' in data_name[-1]: # if a) data_name contains _ and b) genes is the last part of the file name
            dependence = data_name[-1].replace('genes', '').split('.')[-1] # dependence becomes just .filetype 
            data_name = '_'.join(data_name[:-1]) # renames data_name as a string (sep by _) with everything but _genes.filetype

            # full file path
            new_file_name = os.path.join(data_dir, f"{data_name}.{dependence}")
            data_series.append(new_file_name) # add back file type but remove genes

            # create prefix series
            #name_arr = file_name.split() 
            prefix_str = f"{prefix_name}_{data_name}"

            # construct prefix string
            '''if name_arr:
                for ll, part in enumerate(name_arr):
                    if ll == len(name_arr) - 1: # if it's the last(/only) element
                        append_name = part.split('.') # create suffix that is that element split between last word and file type
                        if len(append_name) > 1: # if there is a something.filetype
                            append_name = '_'.join(append_name[:-1]) # recreate suffix into string with just thing without file type
                        else:
                            append_name = part
                    else:
                        append_name = part
                    prefix_str = f"{prefix_str}_{append_name}"
            else:
                prefix_str = f"{prefix_str}_{new_file_name}"'''


            prefix_series.append(prefix_str)

    # ensure unique entries
    data_series = list(set(data_series))
    prefix_series = list(set(prefix_series))

    # create medium_series
    medium_series = [medium] * len(prefix_series)

    return data_series, prefix_series, medium_series
   