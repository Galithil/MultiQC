#!/usr/bin/env python
from __future__ import print_function
import logging
import re
from itertools import chain

from multiqc import config
from multiqc.plots import linegraph, bargraph, scatter, table, heatmap, beeswarm
from multiqc.modules.base_module import BaseMultiqcModule

import numpy as np

from collections import OrderedDict
class slice2OrderedDict(object):
    def __getitem__(self, keys):
        return OrderedDict([(slice.start, slice.stop) for slice in keys])
odict = slice2OrderedDict()

""" MultiQC module to parse output from BBMap """

file_types = {
    'aqhist': {
        'title': 'Average read quality',
        'descr': 'Histogram of average read quality.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'Quality':int,
            'count1':int,
            'fraction1':float,
            'count2':int,
            'fraction2':float
        ],
        'plot_params': {
            'yLog': True,
        }
    },
#basecov
#  - Coverage per base location.
#  - too big to interpret here
    'bhist' : {
        'title': 'Base composition',
        'descr': 'Base composition histogram by position.',
        'help_text': 'Relative composition',
        'cols': odict[
            'Pos':int, 'A':float, 'C':float, 'G':float, 'T':float, 'N':float
        ],
        'plot_params': {}
    },
    'bincov': {
        'title': 'Binned coverage',
        'descr': 'Binned coverage per location (one line per X bases).',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Mean', 'STDev' ],
        'cols': odict[
            'RefName':str,
            'Cov':float,
            'Pos':int,
            'RunningPos':int
        ],
        'plot_params': {},
        'not_implemented': ''
    },
    'bqhist': {
        'title': 'Base quality',
        'descr': 'Quality histogram designed for box plots.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'BaseNum':int,
            'count_1':int, 'min_1':int, 'max_1':int, 'mean_1':float,
            'Q1_1':int, 'med_1':int, 'Q3_1':int, 'LW_1':int, 'RW_1':int,
            'count_2':int, 'min_2':int, 'max_2':int, 'mean_2':float,
            'Q1_2':int, 'med_2':int, 'Q3_2':int, 'LW_2':int, 'RW_2':int
        ],
        'plot_params': {}
    },
    'covhist': {
        'title': 'Coverage histogram',
        'descr': 'Histogram of # occurrences of each depth level.',
        'help_text': 'Placeholder help text.',
        'cols': odict['Coverage':int, 'numBases':int],
        'plot_params': {
            'yLog': True,
        }
    },
    'covstats': {
        'title': 'Coverage stats',
        'descr': 'Per-scaffold coverage info.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'ID':str,
            'Avg_fold':float
        ],
        'extracols': odict[
            'Length':int,
            'Ref_GC':float,
            'Covered_percent':float,
            'Covered_bases':int,
            'Plus_reads':int,
            'Minus_reads':int,
            'Median_fold':int,
            'Read_GC':float,
            'Std_Dev':float
        ],
        'plot_params': {},
        'not_implemented': ''
    },
    'ehist': {
        'title': 'Errors-per-read',
        'descr': 'Errors-per-read histogram.',
        'help_text': 'Placeholder help text.',
        'cols': odict['Errors':int, 'Count':int ],
        'plot_params': {}
    },
    'gchist' : {
        'title': 'GC content',
        'descr': 'Read GC content histogram.',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Mean', 'Median', 'Mode', 'STDev'],
        'cols': odict['GC':float, 'Count':int ],
        'plot_params': {}
    },
    'idhist': {
        'title': 'Identity histogram',
        'descr': 'Histogram of read count versus percent identity.',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Mean_reads', 'Mean_bases',
                   'Median_reads', 'Median_bases',
                   'Mode_reads', 'Mode_bases',
                   'STDev_reads', 'STDev_bases' ],
        'cols': odict[
            'Identity':float, 'Reads':int, 'Bases':int
        ],
        'plot_params': {}
    },
    'ihist': {
        'title': 'Insert sizes',
        'descr': 'Histogram of insert sizes (for paired reads).',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Mean', 'Median', 'STDev', 'PercentOfPairs'],
        'cols': odict['InsertSize':int, 'Count':int ],
        'plot_params': {}
    },
    'indelhist': {
        'title': 'Indel lengths',
        'descr': 'Indel length histogram.',
        'help_text': 'Placeholder help text.',
        'cols': odict['Length':int, 'Deletions':int, 'Insertions':int],
        'plot_params': {}
    },
    'lhist' : {
        'title': 'Read lengths',
        'descr': 'Read length histogram.',
        'help_text': 'Placeholder help text.',
        'cols': odict['Length':int, 'Count':int ],
        'plot_params': {}
    },
    'mhist': {
        'title': 'Match, sub, del, and ins rates',
        'descr': 'Histogram of match, sub, del, and ins rates by read location.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'BaseNum':int,
            'Match1':float, 'Sub1':float, 'Del1':float, 'Ins1':float, 'N1':float, 'Other1':float,
            'Match2':float, 'Sub2':float, 'Del2':float, 'Ins2':float, 'N2':float, 'Other2':float
        ],
        'plot_params': {}
    },
    'qahist': {
        'title': 'Quality accuracy',
        'descr': 'Quality accuracy histogram of error rates versus quality score.',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Deviation', 'DeviationSub' ],
        'cols': odict[
            'Quality':int, 'Match':int, 'Sub':int, 'Ins':int, 'Del':int,
            'TrueQuality':float, 'TrueQualitySub':float
        ],
        'plot_params': {}
    },
    'qhist': {
        'title': 'Quality',
        'descr': 'Quality histogram by position.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'BaseNum':int,
            'Read1_linear':float,
            'Read1_log':float,
            'Read1_measured':float,
            'Read2_linear':float,
            'Read2_log':float,
            'Read2_measured':float
        ],
        'plot_params': {}
    },
    'rpkm': {
        'title': 'RPKM/FPKM',
        'descr': 'Per-scaffold RPKM/FPKM counts.',
        'help_text': 'Placeholder help text.',
        'kvrows': ['File', 'Reads', 'Mapped', 'RefSequences' ],
        'cols': odict[
            'Name':str,
            'Length':int, 'Bases':int, 'Coverage':float,
            'Reads':int, 'RPKM':float, 'Frags':int, 'FPKM':float
        ],
        'plot_params': {},
        'not_implemented': '',
    },
    'statsfile_machine': {
        'title': 'General stats',
        'descr': 'General Stats',
        'help_text': 'Placeholder help text.',
        'cols': [],
        'plot_params': {}
    },
    'statsfile': {
        'title': 'General stats',
        'descr': 'General Stats',
        'help_text': 'Placeholder help text.',
        'cols': [],
        'plot_params': {},
        'not_implemented': ''
    }
}

statsfile_machine_keys = [
    'Reads_Used', 'Bases_Used', 'Reads/sec', 'kBases/sec',
    
    'R1_Mapped_Percent', 'R1_Unambiguous_Percent', 'R1_Mapped_Reads',
    'R1_Unambiguous_Reads', 'Mated_Pairs', 'Bad_Pairs', 'R1_Rescued',
    'Avg_Insert_Size', 'R1_Perfect_Best_Site', 'R1_Semiperfect_Site',
    'R1_Ambiguous_Mapping', 'R1_Low_Quality_Discards', 'R1_Match_Rate',
    'R1_Error_Rate', 'R1_Sub_Rate', 'R1_Del_Rate', 'R1_Ins_Rate',
    'R1_N_Rate', 'R1_Match_Count', 'R1_Error_Count', 'R1_Sub_Count',
    'R1_Del_Count', 'R1_Ins_Count', 'R1_N_Count',

    'R2_Mapped_Percent', 'R2_Unambiguous_Percent', 'R2_Mapped_Reads',
    'R2_Unambiguous_Reads', 'R2_Rescued', 'R2_Perfect_Best_Site',
    'R2_Semiperfect_Site', 'R2_Ambiguous_Mapping', 'R2_Low_Quality_Discards',
    'R2_Match_Rate', 'R2_Error_Rate', 'R2_Sub_Rate', 'R2_Del_Rate',
    'R2_Ins_Rate', 'R2_N_Rate', 'R2_Match_Count', 'R2_Error_Count',
    'R2_Sub_Count', 'R2_Del_Count', 'R2_Ins_Count', 'R2_N_Count',
    'R2_Mapped_Percent',
    
]


# Initialize the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ BBMap module, tries to identify and parse tons of output files
    generated by BBMap
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="BBTools",
            anchor="bbmap",
            href="http://jgi.doe.gov/data-and-tools/bbtools/",
            info="""is a suite of fast multithreaded bioinformatics tools 
            designed for analysis of DNA and RNA sequence data."""
        )

        # Init data dict
        self.mod_data = { key:{} for key in file_types }

        # Find output files
        module_filetypes = [('bbmap/'+ft, ft) for ft in file_types]
        data_found = False
        for module_filetype, file_type in module_filetypes:
            for f in self.find_log_files(module_filetype, filehandles=True):
                if self.parse_logs(file_type, **f):
                    self.add_data_source(f)
                    data_found = True

        if not data_found:
            log.debug("Could not find any data in '%s'", config.analysis_dir)
            raise UserWarning

        for file_type in file_types:
            if len(self.mod_data[file_type]) > 0:
                log.debug("section %s has %d entries", file_type,
                          len(self.mod_data[file_type]))

                self.add_section(
                    name = file_types[file_type]['title']+' ('+file_type+')',
                    anchor =  'bbmap-' + file_type,
                    description = file_types[file_type]['descr'],
                    helptext = file_types[file_type]['help_text'],
                    plot = self.plot_hist(file_type)
                )


    def parse_logs(self, file_type, root, s_name, fn, f, **kw):
        log.debug("Parsing %s/%s", root, fn)
        if not file_type in file_types:
            log.error("Unknown output type '%s'. Error in config?", file_type)
            return False
        log_descr = file_types[file_type]
        if 'not_implemented' in log_descr:
            log.warning("Can't parse '%s' -- implementation missing", file_type)
            return False

        cols = log_descr['cols']
        if isinstance(cols, OrderedDict):
            cols = list(cols.keys())

        kv = {}
        data = {}
        for line_number, line in enumerate(f, start=1):
            line = line.strip().split('\t')
            if line[0][0] == '#':
                # It's a header row
                line[0] = line[0][1:] # remove leading '#'
                    
                if line[0] != cols[0]:
                    # It's not the table header, it must be a key-value row
                    if len(line) !=2:
                        # Not two items? Wrong!
                        log.error("Expected key value pair in %s/%s:%d but found '%s'",
                                  root, s_name, line_number, repr(line))
                        log.error("Table header should begin with '%s'",
                                  cols[0])
                        continue
                    # save key value pair
                    kv[line[0][1:]] = line[1]
                else:
                    # It should be the table header. Verify:
                    if line != cols:
                        if line != cols + list(log_descr['extracols'].keys()):
                            log.error("Table headers do not match those 'on file'. %s != %s",
                                      repr(line), repr(cols))
                        return False
            else:
                if isinstance(log_descr['cols'], OrderedDict):
                    line = [
                        value_type(value)
                        for value_type, value in zip(log_descr['cols'].values(), line)
                    ]
                else:
                    line = list(map(int, line))
                data[line[0]] = line[1:]

        if s_name in self.mod_data[file_type]:
            log.debug("Duplicate sample name found! Overwriting: %s", s_name)

        self.mod_data[file_type][s_name] = {'data':data, 'kv': kv}
        log.debug("Found %s output for sample %s with %d rows",
                  file_type, s_name, len(data))

        return True

    def plot_hist(self, file_type):
        samples = self.mod_data[file_type]

        sumy = sum([int(samples[sample]['data'][x][0])
                    for sample in samples
                    for x in samples[sample]['data']])

        cutoff = sumy * 0.999
        all_x = set()
        for item in sorted(chain(*[samples[sample]['data'].items()
                                   for sample in samples])):
            all_x.add(item[0])
            cutoff -= item[1][0]
            if cutoff < 0:
                xmax = item[0]
                break

        data = {
            sample: {
                x: samples[sample]['data'][x][0] if x in samples[sample]['data'] else 0
                for x in all_x
            }
            for sample in samples
        }

        plot_params = {
                'id': 'bbmap-' + file_type,
                'title':  file_types[file_type]['title'],
                'xmax': xmax
        }
        plot_params.update(file_types[file_type]['plot_params'])
        plot = linegraph.plot(
            data,
            plot_params
        )

        return plot



