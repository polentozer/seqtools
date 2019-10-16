#cSpell: Disable#
GGA_PART_TYPES = {
    'type1': {
        'prefix': 'GCATCGTCTCATCGGAGTCGGTCTCACCCT',
        'suffix': 'AACGAGAGACCAGCAGACCAGAGACGGCAT',
        'info': 'Left side assembly connector'
    },
    'type2': {
        'prefix': 'GCATCGTCTCATCGGTCTCAAACG',
        'suffix': 'TATGAGAGACCTGAGACGGCAT',
        'info': 'Promotor'
    },
    'type3': {
        'prefix': 'GCATCGTCTCATCGGTCTCAT',
        'suffix': 'ATCCAGAGACCTGAGACGGCAT',
        'info': 'CDS'
    },
    'type3a': {
        'prefix': 'GCATCGTCTCATCGGTCTCAT',
        'suffix': 'GGTTCTAGAGACCTGAGACGGCAT',
        'info': 'N-terminal CDS'
    },
    'type3b': {
        'prefix': 'GCATCGTCTCATCGGTCTCATTCT',
        'suffix': 'GGATCCAGAGACCTGAGACGGCAT',
        'info': 'CDS'
    },
    'type3t': {
        'prefix': 'GCATCGTCTCATCGGTCTCAT',
        'suffix': 'GGATCCTGAGACCTGAGACGGCAT',
        'info': 'True type3 CDS (GS linker, no STOP)'
    },
    'type4': {
        'prefix': 'GCATCGTCTCATCGGTCTCAATCC',
        'suffix': 'GCTGAGAGACCTGAGACGGCAT',
        'info': 'Terminator'
    },
    'type4a': {
        'prefix': 'GCATCGTCTCATCGGTCTCAATCC',
        'suffix': 'TGGCAGAGACCTGAGACGGCAT',
        'info': 'C-terminal CDS'
    },
    'type4b': {
        'prefix': 'GCATCGTCTCATCGGTCTCATGGC',
        'suffix': 'GCTGAGAGACCTGAGACGGCAT',
        'info': 'Terminator'
    },
    'type5': {
        'prefix': 'GCATCGTCTCATCGGAGTCGGTCTCAGCTG',
        'suffix': 'TACAAGAGACCAGCAGACCAGAGACGGCAT',
        'info': 'Right side assembly connector'
    },
    'type6': {
        'prefix': 'GCATCGTCTCATCGGTCTCATACA',
        'suffix': 'GAGTAGAGACCTGAGACGGCAT',
        'info': 'Yeast marker'
    },
    'type7': {
        'prefix': 'GCATCGTCTCATCGGTCTCAGAGT',
        'suffix': 'CCGAAGAGACCTGAGACGGCAT',
        'info': "3'-homology or yeast origin"
    },
    'type8': {
        'prefix': 'GCATCGTCTCATCGGTCTCACCGA',
        'suffix': 'CCCTAGAGACCAGAGACGGCAT',
        'info': 'E. coli marker and origin'
    },
    'type8a': {
        'prefix': 'GCATCGTCTCATCGGTCTCACCGA',
        'suffix': 'CAATAGAGACCAGAGACGGCAT',
        'info': 'E. coli marker and origin'
    },
    'type8b': {
        'prefix': 'GCATCGTCTCATCGGTCTCACAAT',
        'suffix': 'CCCTAGAGACCAGAGACGGCAT',
        'info': "5'-homology"
    },
    'typeX': {
        'prefix': 'GCATCGTCTCATCGGTCTCANNNN',
        'suffix': 'NNNNAGAGACCAGAGACGGCAT',
        'info': 'Custom parts'
    }
}

RESTRICTION_ENZYMES = {
    'BsaI': {
        'substrate': 'DNA',
        'recognition': 'GGTCTC',
        'jump': 1,
        'overhang': 4,
        'incubation_temperature': 37,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': True,
            'EcoKI': False
        },
        'description': 'Type 2 restriction enzyme used in modular cloning or MoClo for short. Sticky ends from different BsaI sites may not be compatible. BsaI can be used between 37 and 50 °C.'
    },
    'BsmBI': {
        'substrate': 'DNA',
        'recognition': 'CGTCTC',
        'jump': 1,
        'overhang': 4,
        'incubation_temperature': 55,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': False,
            'EcoKI': False
        },
        'description': 'Type 2 restriction enzyme used in modular cloning or MoClo for short. Sticky ends from different BsmBI sites may not be compatible.'
    },
    'NotI': {
        'substrate': 'DNA',
        'recognition': 'GCGGCCGC',
        'jump': 0,
        'overhang': 4,
        'incubation_temperature': 37,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': False,
            'EcoKI': False
        },
        'description': 'Classic restriction enzyme commonly used in cloning. This particular enzyme has has very rare recognition site.'
    },
    'BpiI': {
        'substrate': 'DNA',
        'recognition': 'GAAGAC',
        'jump': 2,
        'overhang': 4,
        'incubation_temperature': 37,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': False,
            'EcoKI': False
        },
        'description': 'Type 2 restriction enzyme used in modular cloning or MoClo for short. Sticky ends from different BsmBI sites may not be compatible.'
    },
    'EcoRI': {
        'substrate': 'DNA',
        'recognition': 'GAATTC',
        'jump': 0,
        'overhang': 4,
        'incubation_temperature': 37,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': False,
            'EcoKI': False
        },
        'description': 'Classic restriction enzyme commonly used in cloning.'
    }
}

LOGGING_CONFIG = {
    'version': 1,
    'formatters': {
        'simple': {
            'format': '%(asctime)s %(filename)s %(name)s %(levelname)s %(message)s'
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'level': 'DEBUG',
            'formatter': 'simple',
            'stream': 'ext://sys.stdout',
        },
        'file': {
            'class': 'logging.FileHandler',
            'level': 'DEBUG',
            'formatter': 'simple',
            'filename': 'seqtools.log'
        }
    },
    'loggers': {
        '__main__': {
            'level': 'DEBUG',
            'handlers': ['console', 'file'],
            'propagate': False
        }
    },
    'root': {
        'level': 'DEBUG',
        'handlers': ['console']
    }
}

