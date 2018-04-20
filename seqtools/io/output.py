def writer(dictionary, path=None):
    """
    Simple function for writing `.fasta` files.
    Translates dictionary into a fasta file and writes it to a file. If no path is given,
    it returns only translated string.
    """
    string = ""

    for k, v in dictionary.items():
            if v[0]:
                k += "|OPTIMIZED"
                if v[1]:
                    k += "|forced"
                string += "{0}\n{1}\n\n".format(k, v[0])

    if path:
        with open(path, 'w') as file_out:
            file_out.write(string)
    else:
        return string
