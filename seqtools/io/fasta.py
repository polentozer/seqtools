def open_fasta(file_input_paths):
    """
    Function for opening `.fasta` files. Files can contain multiple sequences.
    Returns dictionary of sequences ==> {id: sequence}

    dict = {
        id_1: sequence1,
        id_2: sequence2,
        id_3: sequence3,
        ...
    }

    TODO: change to list:

    list = [
        [id_1, sequence1],
        [id_2, sequence2],
        [id_3, sequence3],
        ...
    ]
    
    """
    sequence_dictionary = {}

    for input_path in file_input_paths:
        with open(input_path, 'r') as input_file:
            data = input_file.readlines()
        
        for line in data:
            temp = line.strip()
            if line[0] == ">":
                sequence_id = temp[:20]
                sequence_dictionary[sequence_id] = ""
            else:
                sequence_dictionary[sequence_id] += temp.upper()

    return sequence_dictionary
