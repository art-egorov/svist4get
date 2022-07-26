from reportlab.pdfgen import canvas
from reportlab.lib.units import cm, mm
from reportlab.lib.colors import black, red
import reportlab.rl_config
import pybedtools

reportlab.rl_config.warnOnMissingFontGlyphs = 0
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import stringWidth
import math
import Bio
from Bio.Seq import Seq
import configs
import argparse
import configs
import os
import sys
from Bio import SeqIO

import svist4get.drawing as drawing
import svist4get.data_processing as data_processing
import svist4get.methods as methods


class Parameters:
    def __init__(self):
        pass

    def parse_command_line(self):
        internal_directory = os.path.join(os.path.dirname(__file__), 'svist4get_data')
        path_to_help = internal_directory + '/help/short_help.txt'
        with open(path_to_help, 'r') as short_help:
            short_help = short_help.read()

        parser = argparse.ArgumentParser(prog='svist4get', description='svist4get commmand line format:',
                                         usage=short_help,
                                         add_help=False)

        parser.add_argument('-tl', dest='transcript_labels', action='append', nargs='*', default=[])

        parser.add_argument('-tls', dest='transcript_label_style',
                            choices=['none', 'name', 'id', 'both', 'empty', 'auto', 'gene_id'], default='empty')

        parser.add_argument("--help", "-h", dest='help_message', action='store_true')

        parser.add_argument("-blp", dest='bedgraph_label_position', choices=['left', 'right', 'center', None],
                            default=None)

        parser.add_argument("--sampledata", "-sampledata", dest='data', action='store_true',
                            help='Creates the \'svist4get\' folder in the current working directory.The folder will contain sample data sets and adjustable configuration file templates.')

        parser.add_argument("--version", "-v", action='version', version='%(prog)s 1.3',
                            help='Show program version')

        parser.add_argument('-t', type=str, metavar='', dest='transcript_id', default=0,
                            help='<id> Visualization of the genomic window of a particular transcript (uses GTF transcript_id)')

        parser.add_argument('-bgc', nargs='*', dest='c_bedgraph_tracks')

        parser.add_argument('-gic', nargs='*', dest='c_regions')

        parser.add_argument('-g', type=str, dest='gene_id', default=0,
                            help='<id>  Visualization of the genomic window of a particular gene (uses GTF gene_id)',
                            metavar='')

        parser.add_argument('-w', dest='window', nargs='*', default=[0],
                            help='tss|tis <upstream> <downstream>  Visualization of a genomic window centered on a transcription start site|translation initiation site of a particular transcript (using GTF transcript_id) with a fixed offsets upstream and downstream.  OR  <contig> <start> <end> Visualization an arbitrary genomic window using genomic coordinates',
                            metavar='')

        parser.add_argument('-fa', dest='fasta_file', type=str, default=0,
                            help='<path/file.fa>  Path to the genome assembly multifasta for the sequence track',
                            metavar='')

        parser.add_argument('-gtf', dest='gtf_file', type=str, default=0,
                            help='<path/file.gtf>   Path to the gtf file with the transcript annotation', metavar='')

        parser.add_argument('-bg', nargs='*', dest='bedgraph', default='', type=str,
                            help='<path/file1.bedGraph>[<path/file2.bedGraph>...]  Paths to bedGraph files for genomic signal tracks.',
                            metavar='')

        parser.add_argument('-c', dest='path_to_config', default='A4_p2', type=str,
                            help='<path/file.cfg> |A4_p1|A4_p2|A4_l   Path to a configuration file or name of the premade config file (can be A4_p1, A4_p2, A4_l) (default = A4_p2)',
                            metavar='')

        parser.add_argument('--debug', '-debug', dest='debug', action='store_true')

        parser.add_argument('--verbose', '-verbose', dest='verbose', default = False, action='store_true')

        parser.add_argument('-o', type=str, dest='output_filename',
                            default='svist4get_output',
                            help='Output files base name. Generates a pair of pdf (vector) and png (raster) files.',
                            metavar='')

        parser.add_argument('-hi', dest='hide_introns', action='store_true',
                            help='Hide intronic segments (default: false)')

        parser.add_argument('-rc', dest='revcomp_transform', action='store_true',
                            help='Reverse-complement transformation of the genomic window (default: off)')

        parser.add_argument('-bg_log', dest='log_scale', action='store_true')

        parser.add_argument('-bgb', dest='bedgraph_bar',
                            choices=['mean', 'max', 'min', 'median', 'none'], metavar='', default='mean',
                            help=' 1|0  Scaling style for bedGraph tracks (default: %(default)s)')

        parser.add_argument('-nts', dest='show_nt_seq_track', default='auto', action='store_true',
                            help='Show reference genomic sequence with single-nucleotide resolution (default: for short regions e.g. < 160 nts)')

        parser.add_argument('-aas', dest='aa_track_style',
                            choices=['codons', 'tics', 'auto'], metavar='', default='auto',
                            help=' largee|small   Scaling style for aminoacid sequence track: large (for long regions) or small (keep single-aminoacid resolution) (default: auto)')

        parser.add_argument('-it', dest='image_title', type=str, default='',
                            help=' <text>   Image title', metavar='')

        parser.add_argument('-xs', dest='gaxis_tics_step', type=int, default=0,
                            help=' <N>  X axis tics step [nt] (default: auto)', metavar='')

        parser.add_argument('-lb', dest='bedgraph_axis_tics', nargs='*', metavar='', default='auto',
                            help=' <N1>[<N2>...]   Y axis tics levels for a bedGraph track, a space-separated list of values, a separate value for each bedGraph track (default: auto)')

        parser.add_argument('-ys', dest='bedgraph_axis_tics_step', metavar='', nargs='*', type=float, default=0,
                            help=' <N1> [<N2>...]    Step between axis tics for a bedgraph track Y-axis, a space-separated list of values, a separated value for each bedgraph track (default:auto)')

        parser.add_argument('-bul', dest='bedgraph_upper_limit', metavar='', nargs='*', default='auto',
                            help=' <N1> <'
                                 'N2>...| max   Upper limit for a bedGraph track Y-axis, can be: 1. a space-separated list of values: a separate axis limit for each bedGraph track; 2.a maximum reachable value across all visible tracks. (default: reachable maximum for each track separately)')

        parser.add_argument('-bll', dest='bedgraph_lower_limit', metavar='', nargs='*', default='auto',
                            help='<N1> <N2>...  Lower limit for a bedGraph track Y-axis, a space-separated list of value, a separate value required for each bedGraph track (default: 0 for each track)')

        parser.add_argument('-bl', dest='bedgraph_label', metavar='', nargs='*', default='',
                            help='<text1> <text2>... Text labels for a bedGraph tracks, space-separated list of values, a separate value required for each bedGraph track (default: names of the bedGraph files)')

        parser.add_argument('-pbg', dest='paired_bedgraph', action='append', nargs='*', default=[])

        parser.add_argument('-gi', dest='genomic_intervals', metavar='', nargs='*', default=[], action='append',
                            type=str,
                            help='<s-e> [<s1-e1>...]  Additional track showing selected genomic intervals as segment, space-separated list of start-end pairs (on genomic coordinates)')

        parser.add_argument('-gil', dest='genomic_intervals_label', metavar='', nargs='*', default='',
                            help='<text1> <text2>...  Labels for genomic intervals, space-separated list of string, a separate value required for each  genomic intervals track')

        parser.add_argument('-st', dest='selected_transcripts', metavar='', nargs='*', default=0,
                            help=' <transcript_id1> <transcript_id2>... A complete list of transcript IDs to shor (other will be hidden, applicable to -w and -g modes)')

        parser.add_argument('-hrf', dest='highlight_reading_frame', metavar='',
                            help=' <N>   Highlight a particular reading frame (+0, +1 or +2)', default=-1)

        parser.add_argument('-hc', dest='highlight_codon', type=str, metavar='', default='',
                            help='<NtNtNt>   Highlight selected codon on the aminoacid sequence track, e.g. -hc ATT')

        parser.add_argument('-hf', dest='highlight_frame', metavar='', action='append', nargs='*', default=[],
                            help='<x1> <x2> <label> Highlight and label a given segment of the genomic window (global coordinates on the selected contig)')

        arg_short_help = 0
        if len(sys.argv[1:]) == 0:
            args = parser.parse_args()
            args = vars(args)
            args['help_message'] = 1
            arg_short_help = 1

        else:
            args = parser.parse_args()
            args = vars(args)

        if args['help_message'] == 1 and arg_short_help == 0:
            internal_directory = os.path.join(os.path.dirname(__file__), 'svist4get_data')
            path_to_help = internal_directory + '/help/help.txt'
            with open(path_to_help, 'r') as help:
                data = help.read()
                print(data)
            sys.exit()
        elif args['help_message'] == 1 and arg_short_help == 1:
            internal_directory = os.path.join(os.path.dirname(__file__), 'svist4get_data')
            path_to_help = internal_directory + '/help/short_help.txt'
            with open(path_to_help, 'r') as short_help:
                short_help = short_help.read()
                print(short_help)
            sys.exit()

        if args['c_bedgraph_tracks'] == None:
            del args['c_bedgraph_tracks']
        if args['c_regions'] == None:
            del args['c_regions']
        if args['bedgraph_label_position'] == None:
            del args['bedgraph_label_position']
        self.cmd = args

    def load_config(self, path_to_config):
        try:
            current_directory = os.getcwd()
            internal_directory = os.path.join(os.path.dirname(__file__), 'svist4get_data')

            if path_to_config == 'A4_p1' or path_to_config == 'A4_p2' or path_to_config == 'A4_l':
                name_of_config = path_to_config + '.cfg'
                path_to_config = internal_directory + '/' + name_of_config

            else:
                path_to_config = path_to_config

            config = configs.load(path_to_config)
            config = config.get_config()

            config = config['root']

            # ---------
            config['additional_title'] = ''
            config['highlight_codon'] = ''
            config['bedgraph_label'] = ''
            config['genomic_intervals'] = []
            config['genomic_intervals_label'] = ''
            config['selected_transcripts'] = 0
            config['bedgraph'] = ''
            config['transcript_id'] = 0
            config['gene_id'] = 0
            config['window'] = [0]
            config['highlight_reading_frame'] = -1
            config['highlight_frame'] = 0
            config['image_title'] = ''
            config['transcript_labels'] = []
            config['log_scale'] = 0
            config['verbose'] = False
            # ------------
            config['path_to_config'] = path_to_config

            self.config = config
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to parse the specified config file. Please check the config file') from error

    def load_palette(self):
        try:
            current_directory = os.getcwd()
            internal_directory = os.path.join(os.path.dirname(__file__), 'svist4get_data')

            path_to_config = self.config['path_to_config']

            path_to_palette = self.config['palette']
            path_to_dir_with_config = os.path.dirname(path_to_config)

            path_to_palette = os.path.join(path_to_dir_with_config, path_to_palette)

            palette = configs.load(path_to_palette)

            palette = palette.get_config()
            palette = palette['root']
            self.palette = palette
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to load palette. Please check the palette path in the config file') from error

    def load_triplet_code(self):
        try:
            current_directory = os.getcwd()
            internal_directory = os.path.join(os.path.dirname(__file__), 'svist4get_data')

            path_to_config = self.config['path_to_config']

            path_to_triplet_code = self.config['triplet_code']
            path_to_dir_with_config = os.path.dirname(path_to_config)

            path_to_triplet_code = os.path.join(path_to_dir_with_config, path_to_triplet_code)

            self.config['triplet_code'] = path_to_triplet_code

            self.path_to_triplet_code = path_to_triplet_code
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to load triplet code. Please check the triplet code path in the config file') from error

    def update_config(self):
        if self.cmd['transcript_label_style'] == 'empty':
            del self.cmd['transcript_label_style']

        try:
            self.config.update(self.cmd)
        except Exception as error:
            raise methods.Svist4getError('Unable to update config') from error

    def load_font(self):
        try:
            current_directory = os.getcwd()
            internal_directory = os.path.join(os.path.dirname(__file__), 'files')
            path_to_mono_font = self.config['mono_font']
            path_to_regular_font = self.config['regular_font']

            path_to_config = self.config['path_to_config']
            path_to_dir_with_config = os.path.dirname(path_to_config)

            path_to_mono_font = os.path.join(path_to_dir_with_config, path_to_mono_font)
            path_to_regular_font = os.path.join(path_to_dir_with_config, path_to_regular_font)

            paths_to_fonts = dict(regular=path_to_regular_font, mono=path_to_mono_font)

            self.paths_to_fonts = paths_to_fonts
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to load font. Please check the font paths in the config file.') from error

    def add_gtf_data(self, gtf_data):
        self.gtf = gtf_data

    def create_pdf(self):
        try:
            pdf = canvas.Canvas(self.config['output_filename'], pagesize=(self.config['page_width'] * cm, 0 * cm))
            self.pdf = pdf

        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create the pdf file. There might be a problem with the \'reportlab\' python package.') from error

    def initialize(self, path_to_config):
        self.load_config(path_to_config)
        self.load_palette()
        self.load_triplet_code()
        self.load_font()
        self.create_pdf()


class Data_preparation:
    def __init__(self, parameters):
        self.parameters = parameters

    def create_dict(self):
        pass


class Loader:
    def __init__(self, preparated_data):
        self.global_parameters = preparated_data['global_parameters']
        self.adjustable_parameters = preparated_data['adjustable_parameters']
        self.config = preparated_data['config']
        self.palette = preparated_data['palette']
        self.paths_to_fonts = preparated_data['paths_to_fonts']

    def extract_data(self):
        pass

    def translation_of_coordinate(self):
        pass

    def hide_introns(self):
        pass

    def revcomp_transform(self):
        pass

    def create_track(self):
        pass


class Vgrid_data_preparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_dict(self):
        pdf = self.parameters.pdf

        global_parameters = dict(pdf=pdf)
        adjustable_parameters = dict()
        config = self.parameters.config
        palette = self.parameters.palette
        paths_to_fonts = self.parameters.paths_to_fonts
        preparated_data = dict(adjustable_parameters=adjustable_parameters, global_parameters=global_parameters,
                               config=config, palette=palette, paths_to_fonts=paths_to_fonts)
        return (preparated_data)


class Vgrid_loader(Loader):
    def __init__(self, preparated_data):
        super().__init__(preparated_data)

    def extract_data(self):
        pass

    def translation_of_coordinate(self):
        pass

    def revcomp_transform(self):
        pass

    def hide_introns(self):
        pass

    def create_track(self):
        pdf = self.global_parameters['pdf']
        input_data = dict(pdf=pdf)

        track = drawing.Vertical_grid(input_data, self.config, self.palette, self.paths_to_fonts)
        tracks = []
        tracks.append(track)
        return (tracks)


class Title_data_preparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_dict(self):
        start = self.parameters.gtf['start']
        end = self.parameters.gtf['end']
        scaffold = self.parameters.gtf['name_of_scaffold']
        strand = '(+)'
        pdf = self.parameters.pdf
        name_of_transcript = ''

        additional_info = self.parameters.config['image_title']

        global_parameters = dict(start=start, end=end, scaffold=scaffold, pdf=pdf, strand=strand,
                                 additional_info=additional_info)
        adjustable_parameters = dict()
        preparated_data = dict(global_parameters=global_parameters, adjustable_parameters=adjustable_parameters,
                               config=self.parameters.config, palette=self.parameters.palette,
                               paths_to_fonts=self.parameters.paths_to_fonts)

        return (preparated_data)


class Title_loader(Loader):
    def __init__(self, preparated_data):
        super().__init__(preparated_data)

    def extract_data(self):
        pass

    def translation_of_coordinate(self):
        pass

    def revcomp_transform(self):
        self.global_parameters['strand'] = '(-)'

    def hide_introns(self):
        pass

    def create_track(self):
        tracks = []
        input_data = dict()
        input_data.update(self.adjustable_parameters)
        input_data.update(self.global_parameters)

        track = drawing.Title(input_data, self.config, self.palette, self.paths_to_fonts)
        tracks.append(track)
        return (tracks)


class Axis_tics_data_praparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_dict(self):

        adjustable_parameters = []
        for i in self.parameters.gtf['for_transcript']:
            dictonary = dict(exons=i['exons'])
            adjustable_parameters.append(dictonary)

        step = self.parameters.config['gaxis_tics_step']
        end = int(self.parameters.gtf['end'])
        start = int(self.parameters.gtf['start'])
        pdf = self.parameters.pdf

        count = start
        coordinates = []
        for i in range(end - start):
            coordinates.append(count)
            count += 1

        global_parameters = dict(pdf=pdf, step=step, start=start, end=end, coordinates=coordinates)

        preparated_data = dict(global_parameters=global_parameters, adjustable_parameters=adjustable_parameters,
                               config=self.parameters.config, palette=self.parameters.palette,
                               paths_to_fonts=self.parameters.paths_to_fonts)
        return (preparated_data)


class Axis_tics_loader(Loader):

    def __init__(self, preparated_data):
        super().__init__(preparated_data)

    def extract_data(self):
        self.coordinates = self.global_parameters['coordinates']
        self.start = self.global_parameters['start']
        self.end = self.global_parameters['end']

    def translation_of_coordinate(self):
        self.start = int(self.global_parameters['start'])
        self.end = int(self.global_parameters['end'])
        start = self.start
        end = self.end
        self.start = int(start) - int(start)
        self.end = int(end) - int(start)

        for k in range(len(self.adjustable_parameters)):

            exons_row = self.adjustable_parameters[k]['exons'].split(';')
            exons = ''
            if exons_row[0] != '':
                for j in exons_row:
                    i = j.split('-')

                    if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                            self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):

                        exons += str(self.global_parameters['start']) + '-' + str(i[1])

                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                            self.global_parameters['end']):
                        exons += str(i[0]) + '-' + str(i[1])

                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[0]) < int(
                            self.global_parameters['end']) and int(i[1]) >= int(self.global_parameters['end']):

                        exons += str(i[0]) + '-' + str(self.global_parameters['end'])

                        exons += ';'
                    elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                            self.global_parameters['end']):
                        exons += str(self.global_parameters['start']) + '-' + (str(self.global_parameters['end']))
                        exons += ';'
                new_exon = ''
                for i in exons.split(';'):
                    if i != '':
                        new_exon += i
                        if i != exons.split(';')[::-1][1]:
                            new_exon += ';'
                exons = new_exon
            else:
                self.exons = ''

            exons = exons.split(';')
            correct_exons = ''

            if len(exons) != 0:
                for j in exons:
                    if j != '':
                        i = j.split('-')
                        s = int(i[0])
                        e = int(i[1])
                        left = s - int(start)
                        right = e - int(start)
                        correct_exons += (str(left) + '-' + str(right))
                        if j != exons[::-1][0]:
                            correct_exons += ';'
                self.exons = correct_exons

            self.adjustable_parameters[k]['exons'] = self.exons.split(';')

    def hide_introns(self):

        line = []
        for i in range(self.end - self.start):
            line.append(i)

        for i in range(len(self.adjustable_parameters)):
            exons = self.adjustable_parameters[i]['exons']
            new_line = []
            for j in exons:
                if j != '':
                    j = j.split('-')
                    start = int(j[0])
                    end = int(j[1])
                    new_line += line[start:end]

            new_line = sorted(new_line)

            introns = []
            count = 0
            if new_line == []:
                introns.append(str(self.start) + '-' + str(self.end))

            for k in range(len(new_line)):
                if count == 0 and new_line[0] != self.start:
                    count += 1
                    introns.append(str(self.start) + '-' + str(new_line[0]))

                if k != len(new_line) - 1:
                    if new_line[k] != (new_line[k + 1] - 1):
                        start = int(new_line[k]) + 1
                        end = new_line[k + 1]
                        introns.append(str(start) + '-' + str(end))
                else:
                    if new_line[len(new_line) - 1] != self.end - 1:
                        start = int(new_line[len(new_line) - 1]) + 1
                        end = self.end
                        introns.append(str(start) + '-' + str(end))
            self.adjustable_parameters[i]['introns'] = introns

        intersect_introns = []
        for i in range(len(self.adjustable_parameters)):
            if i == 0:
                intersect_introns += set(self.adjustable_parameters[i]['introns'])

            if i == 0:

                adjust_param = []
                for k in self.adjustable_parameters[i]['introns']:
                    k = k.split('-')
                    start = int(k[0])
                    end = int(k[1])
                    for l in range(end - start):
                        adjust_param.append(start)
                        start += 1
                intersect_introns1 = adjust_param

            if i != 0:
                adjust_param = []

                for p in self.adjustable_parameters[i]['introns']:
                    p = p.split('-')
                    start = int(p[0])
                    end = int(p[1])
                    for q in range(end - start):
                        adjust_param.append(start)
                        start += 1

            intersect_introns1 = (set(intersect_introns1) & set(adjust_param))

        intersect_introns4 = []
        for o in intersect_introns1:
            intersect_introns4.append(o)
        intersect_introns2 = intersect_introns4

        intersect_introns3 = []
        for z in range(len(intersect_introns2) - 1):
            if z == 0:
                start = int(intersect_introns2[0])
            if intersect_introns2[z] != intersect_introns2[z + 1] - 1:
                end = int(intersect_introns2[z]) + 1
                intersect_introns3.append(str(start) + '-' + str(end))
                start = int(intersect_introns2[z + 1])
            if z == len(intersect_introns2) - 1 and end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1]) + 1

                intersect_introns3.append(str(start) + '-' + str(end))
            if z == (len(intersect_introns2) - 2) and self.end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1] + 1)
                intersect_introns3.append(str(start) + '-' + str(end))

        intersect_introns = intersect_introns3

        true_exon = []
        new_line = []
        for i in intersect_introns:
            if i != '':
                i = i.split('-')
                start = int(i[0])
                end = int(i[1])
                new_line += line[start:end]
        new_line = sorted(new_line)
        count = 0
        if new_line == []:
            true_exon.append(str(self.start) + '-' + str(self.end))
        for i in range(len(new_line)):
            if count == 0 and new_line[0] != self.start:
                count += 1
                true_exon.append(str(self.start) + '-' + str(new_line[0]))
            if i != len(new_line) - 1:
                if new_line[i] != (new_line[i + 1] - 1):
                    start = int(new_line[i]) + 1
                    end = new_line[i + 1]
                    true_exon.append(str(start) + '-' + str(end))
            else:
                if new_line[len(new_line) - 1] != self.end - 1:
                    start = int(new_line[len(new_line) - 1]) + 1
                    end = self.end
                    true_exon.append(str(start) + '-' + str(end))

        new_coordinates = []
        if true_exon != []:

            for i in true_exon:
                SandE = i.split('-')
                start = int(SandE[0])
                end = int(SandE[1])
                new_coordinates += self.coordinates[start:end]
        self.coordinates = new_coordinates

    def revcomp_transform(self):

        self.reverse = 1

        length = len(self.coordinates)

        for k in range(len(self.adjustable_parameters)):

            new_exon = ''
            self.exons = self.adjustable_parameters[k]['exons']
            for i in self.exons:
                if i != '':
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])

                    new_start = (length - end)
                    new_end = (length - start)
                    new_exon += str(new_start) + '-' + str(new_end)
                    if i != self.exons[::-1][0]:
                        new_exon += ';'
            self.adjustable_parameters[k]['exons'] = new_exon.split(';')

        self.coordinates = self.coordinates[::-1]

    def create_track(self):
        tracks = []
        pdf = self.global_parameters['pdf']
        coordinates = self.coordinates
        step = self.global_parameters['step']
        input_data = dict(pdf=pdf, coordinates=coordinates, step=step)

        track = drawing.Axis_tics(input_data, self.config, self.palette, self.paths_to_fonts)
        tracks.append(track)
        return (tracks)


class Nt_seq_data_preparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)
        self.fasta_file = parameters.config['fasta_file']

    def create_dict(self):
        start = self.parameters.gtf['start']
        end = self.parameters.gtf['end']
        pdf = self.parameters.pdf
        scaffold = self.parameters.gtf['name_of_scaffold']

        coordinates = []
        count = start
        for i in range(end - start):
            coordinates.append(count)
            count += 1

        adjustable_parameters = []
        for i in self.parameters.gtf['for_transcript']:
            dictonary = dict(exons=i['exons'], start_codon=i['start_codon'])
            adjustable_parameters.append(dictonary)

        global_parameters = dict(scaffold=scaffold, start=start, end=end, pdf=pdf, coordinates=coordinates,
                                 fasta_file=self.fasta_file)

        preparated_data = dict(global_parameters=global_parameters, adjustable_parameters=adjustable_parameters,
                               config=self.parameters.config, palette=self.parameters.palette,
                               paths_to_fonts=self.parameters.paths_to_fonts)

        return (preparated_data)


class Nt_seq_loader(Loader):
    def __init__(self, preparated_data):
        super().__init__(preparated_data)
        self.path_to_fasta = self.global_parameters['fasta_file']

    def extract_data(self):
        name_of_file = self.path_to_fasta
        name_of_scaffold = self.global_parameters['scaffold']
        start = self.global_parameters['start']
        end = self.global_parameters['end']
        sequnce_obj = data_processing.SequenceExtraction(name_of_file, name_of_scaffold, start, end)
        self.sequence = sequnce_obj.extract().upper()

    def translation_of_coordinate(self):
        for k in range(len(self.adjustable_parameters)):
            self.start = int(self.global_parameters['start'])
            self.end = int(self.global_parameters['end'])
            start = self.start
            end = self.end
            self.start = int(start) - int(start)
            self.end = int(end) - int(start)
            exons_row = self.adjustable_parameters[k]['exons'].split(';')
            exons = ''
            if exons_row[0] != '':
                for j in exons_row:
                    i = j.split('-')

                    if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                            self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):

                        exons += str(self.global_parameters['start']) + '-' + str(i[1])

                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                            self.global_parameters['end']):
                        exons += str(i[0]) + '-' + str(i[1])

                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[0]) < int(
                            self.global_parameters['end']) and int(i[1]) >= int(self.global_parameters['end']):

                        exons += str(i[0]) + '-' + str(self.global_parameters['end'])

                        exons += ';'
                    elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                            self.global_parameters['end']):
                        exons += str(self.global_parameters['start']) + '-' + (str(self.global_parameters['end']))
                        exons += ';'
                new_exon = ''
                for i in exons.split(';'):
                    if i != '':
                        new_exon += i
                        if i != exons.split(';')[::-1][1]:
                            new_exon += ';'
                exons = new_exon
            else:
                self.exons = ''

            exons = exons.split(';')
            correct_exons = ''

            if len(exons) != 0:
                for j in exons:
                    if j != '':
                        i = j.split('-')
                        s = int(i[0])
                        e = int(i[1])
                        left = s - int(start)
                        right = e - int(start)
                        correct_exons += (str(left) + '-' + str(right))
                        if j != exons[::-1][0]:
                            correct_exons += ';'
                self.exons = correct_exons

            self.adjustable_parameters[k]['exons'] = self.exons.split(';')

    def hide_introns(self):

        line = []
        for i in range(self.end - self.start):
            line.append(i)

        for i in range(len(self.adjustable_parameters)):
            exons = self.adjustable_parameters[i]['exons']
            new_line = []
            for j in exons:
                if j != '':
                    j = j.split('-')
                    start = int(j[0])
                    end = int(j[1])
                    new_line += line[start:end]

            new_line = sorted(new_line)

            introns = []
            count = 0
            if new_line == []:
                introns.append(str(self.start) + '-' + str(self.end))

            for k in range(len(new_line)):
                if count == 0 and new_line[0] != self.start:
                    count += 1
                    introns.append(str(self.start) + '-' + str(new_line[0]))

                if k != len(new_line) - 1:
                    if new_line[k] != (new_line[k + 1] - 1):
                        start = int(new_line[k]) + 1
                        end = new_line[k + 1]
                        introns.append(str(start) + '-' + str(end))
                else:
                    if new_line[len(new_line) - 1] != self.end - 1:
                        start = int(new_line[len(new_line) - 1]) + 1
                        end = self.end
                        introns.append(str(start) + '-' + str(end))
            self.adjustable_parameters[i]['introns'] = introns

        intersect_introns = []
        for i in range(len(self.adjustable_parameters)):
            if i == 0:
                intersect_introns += set(self.adjustable_parameters[i]['introns'])

            intersect_introns = (set(intersect_introns) & set(self.adjustable_parameters[i]['introns']))

            if i == 0:

                adjust_param = []
                for k in (self.adjustable_parameters[i]['introns']):
                    k = k.split('-')
                    start = int(k[0])
                    end = int(k[1])
                    for l in range(end - start):
                        adjust_param.append(start)
                        start += 1
                intersect_introns1 = adjust_param

            if i != 0:
                adjust_param = []

                for p in (self.adjustable_parameters[i]['introns']):
                    p = p.split('-')
                    start = int(p[0])
                    end = int(p[1])
                    for q in range(end - start):
                        adjust_param.append(start)
                        start += 1

            intersect_introns1 = (set(intersect_introns1) & set(adjust_param))

        intersect_introns4 = []
        for o in intersect_introns1:
            intersect_introns4.append(o)
        intersect_introns2 = intersect_introns4

        intersect_introns3 = []
        for z in range(len(intersect_introns2) - 1):
            if z == 0:
                start = int(intersect_introns2[0])
            if intersect_introns2[z] != intersect_introns2[z + 1] - 1:
                end = int(intersect_introns2[z]) + 1
                intersect_introns3.append(str(start) + '-' + str(end))
                start = int(intersect_introns2[z + 1])
            if z == len(intersect_introns2) - 1 and end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1]) + 1

                intersect_introns3.append(str(start) + '-' + str(end))
            if z == (len(intersect_introns2) - 2) and self.end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1] + 1)
                intersect_introns3.append(str(start) + '-' + str(end))

        intersect_introns = intersect_introns3

        true_exon = []
        new_line = []
        for i in intersect_introns:
            if i != '':
                i = i.split('-')
                start = int(i[0])
                end = int(i[1])
                new_line += line[start:end]
        new_line = sorted(new_line)
        count = 0
        if new_line == []:
            true_exon.append(str(self.start) + '-' + str(self.end))
        for i in range(len(new_line)):
            if count == 0 and new_line[0] != self.start:
                count += 1
                true_exon.append(str(self.start) + '-' + str(new_line[0]))
            if i != len(new_line) - 1:
                if new_line[i] != (new_line[i + 1] - 1):
                    start = int(new_line[i]) + 1
                    end = new_line[i + 1]
                    true_exon.append(str(start) + '-' + str(end))
            else:
                if new_line[len(new_line) - 1] != self.end - 1:
                    start = int(new_line[len(new_line) - 1]) + 1
                    end = self.end
                    true_exon.append(str(start) + '-' + str(end))

        row_sequence = self.sequence
        sequence = ''
        for i in true_exon:
            SandE = i.split('-')
            start = int(SandE[0])
            end = int(SandE[1])
            sequence += row_sequence[start:end]
        self.sequence = sequence

    def revcomp_transform(self):

        length = len(self.sequence)

        for k in range(len(self.adjustable_parameters)):

            new_exon = ''
            self.exons = self.adjustable_parameters[k]['exons']
            for i in self.exons:
                if i != '':
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])

                    new_start = (length - end)
                    new_end = (length - start)
                    new_exon += str(new_start) + '-' + str(new_end)
                    if i != self.exons[::-1][0]:
                        new_exon += ';'
            self.adjustable_parameters[k]['exons'] = new_exon.split(';')

        self.reverse = 1
        row_sequence = Seq(self.sequence)
        sequence = row_sequence.reverse_complement()
        self.sequence = str(sequence)

    def create_track(self):
        tracks = []
        sequence = self.sequence
        pdf = self.global_parameters['pdf']
        input_data = dict(sequence=sequence, pdf=pdf)

        track = drawing.Nucleotide_sequence(input_data, self.config, self.palette, self.paths_to_fonts)
        tracks.append(track)

        return (tracks)


class Transcript_struct_data_preparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_dict(self):
        start = self.parameters.gtf['start']
        end = self.parameters.gtf['end']
        pdf = self.parameters.pdf
        scaffold = self.parameters.gtf['name_of_scaffold']
        transcript_labels = self.parameters.config['transcript_labels']
        coordinates = []
        count = start
        for i in range(end - start):
            coordinates.append(count)
            count += 1
        adjustable_parameters = []

        for i in self.parameters.gtf['for_transcript']:
            dictonary = dict(gene_id = i['gene_id'] ,coordinates=coordinates, exons=i['exons'], CDS=i['CDS'], strand=i['strand'],
                             name=i['name_of_transcript'], start=i['start_of_transcript'], end=i['end_of_transcript'],
                             transcript_name=i['transcript_name'])
            adjustable_parameters.append(dictonary)

        global_parameters = dict(scaffold=scaffold, start=start, end=end, pdf=pdf, transcript_labels=transcript_labels)

        preparated_data = dict(global_parameters=global_parameters, adjustable_parameters=adjustable_parameters,
                               config=self.parameters.config, palette=self.parameters.palette,
                               paths_to_fonts=self.parameters.paths_to_fonts)

        return (preparated_data)


class Transcript_struct_loader(Loader):
    def __init__(self, preparated_data):
        super().__init__(preparated_data)

    def extract_data(self):

        name_of_scaffold = self.global_parameters['scaffold']
        start = self.global_parameters['start']
        end = self.global_parameters['end']
        self.reverse = 0

        for i in range(len(self.adjustable_parameters)):
            self.start_t = self.adjustable_parameters[i]['start']
            self.end_t = self.adjustable_parameters[i]['end']
            self.transcript = str(self.start_t) + '-' + str(self.end_t)
            self.adjustable_parameters[i]['transcript'] = self.transcript
            self.coordinates = self.adjustable_parameters[0]['coordinates']

    def translation_of_coordinate(self):

        for k in range(len(self.adjustable_parameters)):
            self.start = int(self.global_parameters['start'])
            self.end = int(self.global_parameters['end'])
            start = self.start
            end = self.end
            self.start = int(start) - int(start)
            self.end = int(end) - int(start)
            exons_row = self.adjustable_parameters[k]['exons'].split(';')
            exons = ''
            if exons_row[0] != '':
                for j in exons_row:
                    i = j.split('-')
                    if int(i[0]) > int(i[1]):
                        a = i[1]

                        i[1] = i[0]
                        i[0] = a

                    if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                            self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):

                        exons += str(self.global_parameters['start']) + '-' + str(i[1])

                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                            self.global_parameters['end']):
                        exons += str(i[0]) + '-' + str(i[1])

                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[0]) < int(
                            self.global_parameters['end']) and int(i[1]) >= int(self.global_parameters['end']):

                        exons += str(i[0]) + '-' + str(self.global_parameters['end'])

                        exons += ';'
                    elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                            self.global_parameters['end']):
                        exons += str(self.global_parameters['start']) + '-' + (str(self.global_parameters['end']))
                        exons += ';'
                new_exon = ''
                for i in exons.split(';'):
                    if i != '':
                        new_exon += i
                        if i != exons.split(';')[::-1][1]:
                            new_exon += ';'
                exons = new_exon
            else:
                self.exons = ''
            exons = exons.split(';')
            correct_exons = ''
            if len(exons) != 0:
                for j in exons:
                    if j != '':

                        i = j.split('-')
                        if int(i[0]) > int(i[1]):
                            a = i[1]
                            i[1] = i[0]
                            i[0] = a
                        s = int(i[0])
                        e = int(i[1])
                        left = s - int(start)
                        right = e - int(start)
                        correct_exons += (str(left) + '-' + str(right))
                        if j != exons[::-1][0]:
                            correct_exons += ';'
                self.exons = correct_exons
            self.adjustable_parameters[k]['exons'] = self.exons.split(';')

            CDS_row = self.adjustable_parameters[k]['CDS'].split(';')
            CDS = ''
            if CDS_row[0] != '':
                for j in CDS_row:
                    i = j.split('-')
                    if int(i[0]) > int(i[1]):
                        a = i[1]
                        i[1] = i[0]
                        i[0] = a
                    if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                            self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):

                        CDS += str(self.global_parameters['start']) + '-' + str(i[1])

                        CDS += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                            self.global_parameters['end']):
                        CDS += str(i[0]) + '-' + str(i[1])

                        CDS += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[0]) < int(
                            self.global_parameters['end']) and int(i[1]) >= int(self.global_parameters['end']):

                        CDS += str(i[0]) + '-' + str(self.global_parameters['end'])

                        CDS += ';'
                    elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                            self.global_parameters['end']):
                        CDS += str(self.global_parameters['start']) + '-' + (str(self.global_parameters['end']))
                        CDS += ';'

                new_CDS = ''
                for i in CDS.split(';'):
                    if i != '':
                        new_CDS += i
                        if i != CDS.split(';')[::-1][1]:
                            new_CDS += ';'
                CDS = new_CDS

            else:
                CDS = ''
            CDS = CDS.split(';')
            correct_CDS = ''

            if len(CDS[0]) != 0:
                for j in CDS:
                    i = j.split('-')
                    s = int(i[0])
                    e = int(i[1])
                    left = s - int(start)
                    right = e - int(start)
                    correct_CDS += (str(left) + '-' + str(right))
                    if j != CDS[::-1][0]:
                        correct_CDS += ';'
                self.CDS = correct_CDS
            else:
                self.CDS = ''

            self.adjustable_parameters[k]['CDS'] = self.CDS.split(';')
            transcript_row = self.adjustable_parameters[k]['transcript'].split(';')
            if transcript_row != ['0-0']:
                transcript = ''
                if transcript_row != '':

                    for j in transcript_row:
                        i = j.split('-')
                        if int(i[0]) > int(i[1]):
                            a = i[1]
                            i[1] = i[0]
                            i[0] = a

                        if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                                self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):
                            transcript += str(self.global_parameters['start']) + '-' + str(i[1])
                            # transcript += ';'
                        elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                                self.global_parameters['end']):
                            transcript += str(i[0]) + '-' + str(i[1])
                            # transcript += ';'
                        elif int(i[0]) > int(self.global_parameters['start']) and int(i[0]) < int(
                                self.global_parameters['end']) and int(i[1]) > int(self.global_parameters['end']):
                            transcript += str(i[0]) + '-' + str(self.global_parameters['end'])
                            # transcript += ';'
                        elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                                self.global_parameters['end']):
                            transcript += str(self.global_parameters['start']) + '-' + (
                                str(self.global_parameters['end']))
                            # transcript += ';'
                    # new_transcript = ''
                    # for i in transcript.split(';'):
                    # if i!= '':
                    # new_transcript += i
                    # if i!=transcript.split(';')[::-1][1]:
                    # new_transcript+=';'
                    # transcript = new_transcript
                transcript = transcript.split(';')
                correct_transcript = ''
                for j in transcript:
                    if j != '':
                        i = j.split('-')
                        s = int(i[0])
                        e = int(i[1])
                        left = s - int(start)
                        right = e - int(start)
                        correct_transcript += (str(left) + '-' + str(right))
                        # if j!= transcript[::-1][0]:
                        # correct_transcript += ';'
                self.transcript = correct_transcript
                self.adjustable_parameters[k]['transcript'] = self.transcript
            else:
                self.adjustable_parameters[k]['transcript'] == ['0-0']

    def hide_introns(self):

        line = []
        for i in range(self.end - self.start):
            line.append(i)

        for i in range(len(self.adjustable_parameters)):
            exons = self.adjustable_parameters[i]['exons']
            new_line = []
            for j in exons:
                if j != '':
                    j = j.split('-')
                    start = int(j[0])
                    end = int(j[1])
                    new_line += line[start:end]

            new_line = sorted(new_line)
            introns = []
            count = 0
            if new_line == []:
                introns.append(str(self.start) + '-' + str(self.end))

            for k in range(len(new_line)):
                if count == 0 and new_line[0] != self.start:
                    count += 1
                    introns.append(str(self.start) + '-' + str(new_line[0]))

                if k != len(new_line) - 1:
                    if new_line[k] != (new_line[k + 1] - 1):
                        start = int(new_line[k]) + 1
                        end = new_line[k + 1]
                        introns.append(str(start) + '-' + str(end))
                else:
                    if new_line[len(new_line) - 1] != self.end - 1:
                        start = int(new_line[len(new_line) - 1]) + 1
                        end = self.end
                        introns.append(str(start) + '-' + str(end))
            self.adjustable_parameters[i]['introns'] = introns
        intersect_introns = []
        for i in range(len(self.adjustable_parameters)):
            if i == 0:
                intersect_introns += set(self.adjustable_parameters[i]['introns'])

            intersect_introns = (set(intersect_introns) & set(self.adjustable_parameters[i]['introns']))

            if i == 0:

                adjust_param = []
                for k in self.adjustable_parameters[i]['introns']:
                    k = k.split('-')
                    start = int(k[0])
                    end = int(k[1])
                    for l in range(end - start):
                        adjust_param.append(start)
                        start += 1
                intersect_introns1 = adjust_param

            if i != 0:
                adjust_param = []

                for p in self.adjustable_parameters[i]['introns']:
                    p = p.split('-')
                    start = int(p[0])
                    end = int(p[1])
                    for q in range(end - start):
                        adjust_param.append(start)
                        start += 1

            intersect_introns1 = (set(intersect_introns1) & set(adjust_param))

        intersect_introns4 = []
        for o in intersect_introns1:
            intersect_introns4.append(o)
        intersect_introns2 = intersect_introns4

        intersect_introns3 = []
        for z in range(len(intersect_introns2) - 1):
            if z == 0:
                start = int(intersect_introns2[0])
            if intersect_introns2[z] != intersect_introns2[z + 1] - 1:
                end = int(intersect_introns2[z]) + 1
                intersect_introns3.append(str(start) + '-' + str(end))
                start = int(intersect_introns2[z + 1])
            if z == len(intersect_introns2) - 1 and end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1]) + 1

                intersect_introns3.append(str(start) + '-' + str(end))
            if z == (len(intersect_introns2) - 2) and self.end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1] + 1)
                intersect_introns3.append(str(start) + '-' + str(end))

        intersect_introns = intersect_introns3

        true_exon = []
        new_line = []
        for i in intersect_introns:
            if i != '':
                i = i.split('-')
                start = int(i[0])
                end = int(i[1])
                new_line += line[start:end]
        new_line = sorted(new_line)
        count = 0
        if new_line == []:
            true_exon.append(str(self.start) + '-' + str(self.end))
        for i in range(len(new_line)):
            if count == 0 and new_line[0] != self.start:
                count += 1
                true_exon.append(str(self.start) + '-' + str(new_line[0]))
            if i != len(new_line) - 1:
                if new_line[i] != (new_line[i + 1] - 1):
                    start = int(new_line[i]) + 1
                    end = new_line[i + 1]
                    true_exon.append(str(start) + '-' + str(end))
            else:
                if new_line[len(new_line) - 1] != self.end - 1:
                    start = int(new_line[len(new_line) - 1]) + 1
                    end = self.end
                    true_exon.append(str(start) + '-' + str(end))

        self.exons = true_exon
        for k in range(len(self.adjustable_parameters)):
            row_coordinates = self.coordinates
            coordinates = []

            for i in self.exons:
                SandE = i.split('-')
                start = int(SandE[0])
                end = int(SandE[1])
                coordinates += row_coordinates[start:end]
            self.adjustable_parameters[k]['coordinates'] = coordinates

            exist_exon = []

            for i in range(self.end - self.start):
                exist_exon.append(1)

            for i in self.exons:
                SandE = i.split('-')
                start = int(SandE[0])
                end = int(SandE[1])
                start_for_end = []
                count = start

                for j in range((end - start)):
                    start_for_end.append(count)
                    count += 1

                for j in start_for_end:
                    exist_exon[j] = 0

            count = 1
            for i in range(len(exist_exon)):
                if exist_exon[i] == 1:
                    exist_exon[i] = count
                    count += 1

            for i in range(len(exist_exon)):
                if i != 0:
                    if exist_exon[i] == 0:
                        exist_exon[i] = exist_exon[i - 1]
            exist_exon = exist_exon[::-1]
            exist_exon.append(0)
            exist_exon = exist_exon[::-1]

            exons = self.adjustable_parameters[k]['exons']

            if exons != ['']:
                new_exons = ''
                for i in exons:
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])
                    new_start = start - exist_exon[start]
                    new_end = end - exist_exon[end]
                    new_exons += str(new_start) + '-' + str(new_end)
                    if i != exons[::-1][0]:
                        new_exons += ';'

                self.adjustable_parameters[k]['exons'] = new_exons.split(';')

            CDS = self.adjustable_parameters[k]['CDS']
            new_CDS = ''
            if CDS[0] != '':
                for i in CDS:
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])
                    new_start = start - exist_exon[start]
                    new_end = end - exist_exon[end]
                    new_CDS += str(new_start) + '-' + str(new_end)
                    if i != CDS[::-1][0]:
                        new_CDS += ';'

                self.CDS = new_CDS
            else:
                self.CDS = ''
            self.adjustable_parameters[k]['CDS'] = self.CDS.split(';')

            transcript = self.adjustable_parameters[k]['transcript']

            new_transcript = ''
            i = transcript
            if i != '':
                SandE = i.split('-')
                start = int(SandE[0])
                end = int(SandE[1])
                new_start = start - exist_exon[start]
                new_end = end - exist_exon[end]
                new_transcript += str(new_start) + '-' + str(new_end)

            self.adjustable_parameters[k]['transcript'] = new_transcript

    def revcomp_transform(self):
        self.adjustable_parameters = self.adjustable_parameters[::-1]
        for k in range(len(self.adjustable_parameters)):
            self.reverse = 1
            self.adjustable_parameters[k]['coordinates'] = self.adjustable_parameters[k]['coordinates'][::-1]

            length = len(self.adjustable_parameters[k]['coordinates'])
            new_exon = ''
            self.exons = self.adjustable_parameters[k]['exons']
            for i in self.exons:
                if i != '':
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])

                    new_start = (length - end)
                    new_end = (length - start)
                    new_exon += str(new_start) + '-' + str(new_end)
                    if i != self.exons[::-1][0]:
                        new_exon += ';'
            self.adjustable_parameters[k]['exons'] = new_exon.split(';')[::-1]

            new_CDS = ''
            self.CDS = self.adjustable_parameters[k]['CDS']
            if self.CDS[0] != '':
                for i in self.CDS:
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])

                    new_start = (length - end)
                    new_end = (length - start)
                    new_CDS += str(new_start) + '-' + str(new_end)
                    if i != self.CDS[::-1][0]:
                        new_CDS += ';'
            self.adjustable_parameters[k]['CDS'] = new_CDS.split(';')[::-1]

            new_transcript = ''

            self.transcript = self.adjustable_parameters[k]['transcript']
            i = self.transcript
            if i != '':
                SandE = i.split('-')
                start = int(SandE[0])
                end = int(SandE[1])

                new_start = (length - end)
                new_end = (length - start)
                new_transcript += str(new_start) + '-' + str(new_end)

            self.adjustable_parameters[k]['transcript'] = new_transcript
        self.adjustable_parameters = self.adjustable_parameters

    def create_track(self):
        tracks = []
        length = len(self.coordinates)
        self.length_of_sequence = length
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))
        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))

        self.available_width = self.config['page_width'] * cm
        self.width_of_space = ((self.available_width - 2 * self.config['margin'] * cm) / self.length_of_sequence) * 0.1

        width_of_space = self.width_of_space

        self.width_of_box = ((self.available_width - 2 * self.config[
            'margin'] * cm + self.width_of_space) - self.width_of_space * self.length_of_sequence) / self.length_of_sequence
        width_of_box = self.width_of_box

        for k in range(len(self.adjustable_parameters)):
            if k == 0:
                self.adjustable_parameters[k]['line'] = 0



            else:
                lines = []
                for j in range(len(self.adjustable_parameters)):
                    if j < k:
                        regions1 = self.adjustable_parameters[j]['transcript']
                        regions2 = self.adjustable_parameters[k]['transcript']
                        regions1 = regions1.split('-')
                        regions2 = regions2.split('-')
                        start1 = min(int(regions1[0]), int(regions1[1]))
                        end1 = max(int(regions1[0]), int(regions1[1]))
                        start2 = min(int(regions2[0]), int(regions2[1]))
                        end2 = max(int(regions2[0]), int(regions2[1]))

                        self.name1 = self.adjustable_parameters[j]['name']
                        self.name2 = self.adjustable_parameters[k]['name']

                        find1 = 0
                        find2 = 0
                        for i in self.config['transcript_labels']:
                            if self.name1 in i:
                                self.name1 = i[1]
                                find1 = 1
                            if self.name2 in i:
                                self.name2 = i[1]
                                find2 = 1
                        if find1 == 0 and self.config['transcript_label_style'] == 'none':
                            self.name1 = ''

                        if find2 == 0 and self.config['transcript_label_style'] == 'none':
                            self.name2 = ''

                        if self.config['transcript_label_style'] == 'name':
                            self.name1 = self.adjustable_parameters[j]['transcript_name']
                            self.name2 = self.adjustable_parameters[k]['transcript_name']
                            if self.name1 == '' or self.name1 == ' ':
                                self.name1 = 'N/A'
                            if self.name2 == '' or self.name2 == ' ':
                                self.name2 = 'N/A'

                        if self.config['transcript_label_style'] == 'gene_id':
                            self.name1 = self.adjustable_parameters[j]['gene_id']
                            self.name2 = self.adjustable_parameters[k]['gene_id']
                            if self.name1 == '' or self.name1 == ' ':
                                self.name1 = 'N/A'
                            if self.name2 == '' or self.name2 == ' ':
                                self.name2 = 'N/A'

                        if self.config['transcript_label_style'] == 'both':
                            if self.adjustable_parameters[j]['transcript_name'] != '' and self.name1 not in \
                                    self.adjustable_parameters[j]['transcript_name']:
                                self.name1 = self.name1 + ' (' + self.adjustable_parameters[j]['transcript_name'] + ')'

                            if self.adjustable_parameters[k]['transcript_name'] != '' and self.name2 not in \
                                    self.adjustable_parameters[k]['transcript_name']:
                                self.name2 = self.name2 + ' (' + self.adjustable_parameters[k]['transcript_name'] + ')'

                        if self.config['transcript_label_style'] == 'auto':
                            if self.adjustable_parameters[j]['transcript_name'] != '':
                                self.name1 = self.name1 + ' (' + self.adjustable_parameters[j]['transcript_name'] + ')'

                            if self.adjustable_parameters[k]['transcript_name'] != '':
                                self.name2 = self.name2 + ' (' + self.adjustable_parameters[k]['transcript_name'] + ')'

                        if self.config['transcript_label_style'] == 'none':
                            end_string1 = 0
                            start_string2 = 1

                        self.global_parameters['pdf'].setFont('regular', self.config['transcript_id_font_size'])
                        width_of_string1 = stringWidth(self.name1, 'regular', self.config['transcript_id_font_size'])
                        width_of_string2 = stringWidth(self.name2, 'regular', self.config['transcript_id_font_size'])
                        width1 = (end1 - start1) * width_of_box
                        width2 = (end2 - start2) * width_of_box
                        first = []

                        center1 = ((end1 + start1) * width_of_box) / 2
                        center2 = ((end2 + start2) * width_of_box) / 2

                        end_string1 = center1 + width_of_string1 / 2
                        start_string2 = center2 - width_of_string2 / 2

                        if start_string2 < 0:
                            start_string2 = self.config['margin'] * cm
                        if end_string1 < width_of_string1 + self.config['margin'] * cm:
                            end_string1 = width_of_string1 + self.config['margin'] * cm
                        if end_string1 > self.config['page_width'] * cm - self.config['margin'] * cm:
                            end_string1 = self.config['page_width'] * cm - self.config['margin'] * cm
                        if start_string2 + width_of_string2 > (
                                self.config['page_width'] * cm - 3 * self.config['margin'] * cm):
                            start_string2 = self.config['page_width'] * cm - 5 * self.config[
                                'margin'] * cm - width_of_string2

                        l = start1
                        for i in range(end1 - start1):
                            first.append(l)
                            l += 1
                        second = []
                        p = start2
                        for i in range(end2 - start2):
                            second.append(p)
                            p += 1
                        result = list(set(first) & set(second))
                        a1 = min(first)
                        a2 = max(first)
                        b1 = min(second)
                        b2 = max(second)
                        length1 = end1 - start1
                        length2 = end2 - start2
                        dist = min(abs(b1 - a1), abs(b1 - a2), abs(b2 - a1), abs(b2 - a2))

                        plus1 = 0
                        plus2 = 0
                        if width_of_string1 > width1:
                            plus1 = width_of_string1 - width1
                        if width_of_string2 > width2:
                            plus2 = width_of_string2 - width2
                        dist1 = dist * width_of_box

                        if (len(result) == 0) and (end_string1 < start_string2):
                            if self.adjustable_parameters[j]['line'] not in lines:
                                lines.append(self.adjustable_parameters[j]['line'])
                        else:
                            if self.adjustable_parameters[j]['line'] in lines:
                                lines.remove(self.adjustable_parameters[j]['line'])

            if k != 0:
                if lines != []:
                    self.adjustable_parameters[k]['line'] = min(lines)
                else:
                    self.adjustable_parameters[k]['line'] = -1

            if self.adjustable_parameters[k]['line'] == -1:
                allline = []

                for p in range(len(self.adjustable_parameters)):
                    if p < k:
                        allline.append(self.adjustable_parameters[p]['line'])

                self.adjustable_parameters[k]['line'] = max(allline) + 1

        self.new_adjustable_parameters = []
        for p in range(len(self.adjustable_parameters)):
            for j in range(len(self.adjustable_parameters)):
                if self.adjustable_parameters[j]['line'] == p:
                    if self.adjustable_parameters[j] not in self.new_adjustable_parameters:
                        self.new_adjustable_parameters.append(self.adjustable_parameters[j])

        self.adjustable_parameters = self.new_adjustable_parameters

        num_line = []
        for p in range(len(self.adjustable_parameters)):
            num_line.append(self.adjustable_parameters[p]['line'])

        max_line = max(num_line)

        for k in range(len(self.adjustable_parameters)):
            max_line = max_line
            exons = self.adjustable_parameters[k]['exons']
            coordinates = self.adjustable_parameters[k]['coordinates']
            CDS = self.adjustable_parameters[k]['CDS']
            # end = len(sequence) + 1
            # start = 0
            i = self.adjustable_parameters[k]['transcript']
            if i != '':
                i = i.split('-')
                left = i[0]
                right = i[1]
            else:
                left = len(coordinates)
                right = 0
            line = self.adjustable_parameters[k]['line']
            new_line = 0

            if k != len(self.adjustable_parameters) - 1:
                if line != self.adjustable_parameters[k + 1]['line']:
                    new_line = 1
            else:
                new_line = 1

            delete = str(0) + '-' + str(left) + ';' + str(right) + '-' + str(len(coordinates))
            delete = delete.split(';')
            i = self.adjustable_parameters[k]['transcript'].split('-')
            pdf = self.global_parameters['pdf']
            reverse = self.reverse
            '''
            if self.config['transcript_label'] == 'id':
                name = self.adjustable_parameters[k]['name']
            elif self.config['transcript_label'] == 'name':
                name  = self.adjustable_parameters[k]['transcript_name']
            elif self.config['transcript_label'] == 'both':
                if self.adjustable_parameters[k]['transcript_name'] != '':
                    name = self.adjustable_parameters[k]['transcript_name'] + ' (' + self.adjustable_parameters[k]['name'] + ')'
                else:
                    name = self.adjustable_parameters[k]['name']

            elif self.config['transcript_label'] == 'none':
                name = ''
            '''
            find = 0
            name = self.adjustable_parameters[k]['name']
            for p in self.config['transcript_labels']:
                if name in p:
                    name = p[1]
                    find = 1
            if find == 0 and self.config['transcript_label_style'] == 'none':
                name = ''

            if self.config['transcript_label_style'] == 'name':
                name = self.adjustable_parameters[k]['transcript_name']
                if name == '' or name == ' ':
                    name = 'N/A'
            if self.config['transcript_label_style'] == 'gene_id':
                name = self.adjustable_parameters[k]['gene_id']
                if name == '' or name == ' ':
                    name = 'N/A'

            if self.config['transcript_label_style'] == 'auto':
                if self.adjustable_parameters[k]['transcript_name'] != '' and name not in \
                        self.adjustable_parameters[k]['transcript_name']:
                    name = name + ' (' + self.adjustable_parameters[k]['transcript_name'] + ')'

            if self.config['transcript_label_style'] == 'both':
                if self.adjustable_parameters[k]['transcript_name'] != '':
                    name = name + ' (' + self.adjustable_parameters[k]['transcript_name'] + ')'

            transcript_labels = self.global_parameters['transcript_labels']
            orientation = self.adjustable_parameters[k]['strand']
            input_data = dict(exons=exons, transcript_labels=transcript_labels, max_line=max_line, CDS=CDS,
                              coordinates=coordinates, new_line=new_line, transcript=i,
                              delete=delete,
                              reverse=reverse, name_of_transcript=name, orientation=orientation, pdf=pdf)
            track = drawing.Transcript_structure(input_data, self.config, self.palette, self.paths_to_fonts)
            tracks.append(track)
        return (tracks)


'''

        
        for k in range(len(self.adjustable_parameters)):
                
            self.adjustable_parameters[k]['intersect'] = 1
            if k != (len(self.adjustable_parameters) - 1):
                regions1 = self.adjustable_parameters[k]['transcript']
                regions2 = self.adjustable_parameters[k + 1]['transcript']
                regions1 = regions1.split('-')
                regions2 = regions2.split('-')
                start1 = min(int(regions1[0]), int(regions1[1]))
                end1 = max(int(regions1[0]), int(regions1[1]))
                start2 = min(int(regions2[0]), int(regions2[1]))
                end2 = max(int(regions2[0]), int(regions2[1]))
                self.name1 = self.adjustable_parameters[k]['name']
                self.name2 = self.adjustable_parameters[k + 1]['name']
                self.global_parameters['pdf'].setFont('regular', self.config['transcript_id_font_size'])
                width_of_string1 = stringWidth(self.name1, 'regular', self.config['transcript_id_font_size'])
                width_of_string2 = stringWidth(self.name1, 'regular', self.config['transcript_id_font_size'])
                width1 = (end1 - start1) * width_of_box
                width2 = (end2 - start2) * width_of_box
                first = []

                l = start1
                for i in range(end1 - start1):
                    first.append(l)
                    l += 1
                second = []
                p = start2
                for i in range(end2 - start2):
                    second.append(p)
                    p += 1
                result = list(set(first) & set(second))
                a1 = min(first)
                a2 = max(first)
                b1 = min(second)
                b2 = max(second)
                length1 = end1 - start1
                length2 = end2 - start2
                dist = min(abs(b1 - a1), abs(b1 - a2), abs(b2 - a1), abs(b2 - a2))

                plus1 = 0
                plus2 = 0
                if width_of_string1 > width1:
                    plus1 = width_of_string1 - width1
                if width_of_string2 > width2:
                    plus2 = width_of_string2 - width2
                dist1 = dist * width_of_box

                if len(result) == 0 and abs(dist) > 0.001 * length and (plus1 + plus2) < dist1:
                    self.adjustable_parameters[k]['intersect'] = 0

        for k in range(len(self.adjustable_parameters)):
            exons = self.adjustable_parameters[k]['exons']
            coordinates = self.adjustable_parameters[k]['coordinates']
            CDS = self.adjustable_parameters[k]['CDS']
            # end = len(sequence) + 1
            # start = 0
            i = self.adjustable_parameters[k]['transcript']
            if i != '':
                i = i.split('-')
                left = i[0]
                right = i[1]
            else:
                left = len(coordinates)
                right = 0
            intersect = self.adjustable_parameters[k]['intersect']
            delete = str(0) + '-' + str(left) + ';' + str(right) + '-' + str(len(coordinates))
            delete = delete.split(';')
            i = self.adjustable_parameters[k]['transcript'].split('-')
            pdf = self.global_parameters['pdf']
            reverse = self.reverse
            name = self.adjustable_parameters[k]['name']
            if name == 'MATa2(fs)':
                intersect = 'first'
            orientation = self.adjustable_parameters[k]['strand']
            input_data = dict(exons=exons, CDS=CDS, coordinates=coordinates, intersect=intersect, transcript=i,
                              delete=delete,
                              reverse=reverse, name_of_transcript=name, orientation=orientation, pdf=pdf)
            track = drawing.Transcript_structure(input_data, self.config, self.palette, self.paths_to_fonts)
            tracks.append(track)
        return (tracks)
'''


class Bedgraph_data_preparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_dict(self):
        start = self.parameters.gtf['start']
        end = self.parameters.gtf['end']
        step = self.parameters.config['bedgraph_axis_tics_step']
        lvl_of_tics = self.parameters.config['bedgraph_axis_tics']
        self.maxmax = 0
        if self.parameters.config['bedgraph_upper_limit'] != 'd':
            if self.parameters.config['bedgraph_upper_limit'][0] == 'max':
                self.parameters.config['bedgraph_upper_limit'] = 'd'
                self.maxmax = 1

        upper_limit = self.parameters.config['bedgraph_upper_limit']
        scaffold = self.parameters.gtf['name_of_scaffold']
        all_exons = []
        for l in range(len(self.parameters.gtf['for_transcript'])):
            exons = self.parameters.gtf['for_transcript'][l]['exons']
            dictonary = dict(exons=exons)
            all_exons.append(dictonary)
        pdf = self.parameters.pdf

        adjustable_parameters = []
        if self.parameters.config['bedgraph_label'] == '':
            self.parameters.config['bedgraph_label'] = []
            for q in self.parameters.config['bedgraph']:
                q = os.path.abspath(q)
                q = os.path.basename(q)
                q = os.path.splitext(q)[0]
                self.parameters.config['bedgraph_label'].append(q)

        for i in range(len(self.parameters.config['bedgraph'])):

            if self.parameters.config['bedgraph_upper_limit'] == 'auto' and self.maxmax == 0:
                upper_limit = 'd'
            elif self.maxmax == 0:
                upper_limit = float(self.parameters.config['bedgraph_upper_limit'][i])
            self.minmin = 0
            if self.parameters.config['bedgraph_lower_limit'] != 'auto':
                if self.parameters.config['bedgraph_lower_limit'][0] == 'min':
                    self.parameters.config['bedgraph_lower_limit'] = 'auto'
                    self.minmin = 1

            if self.parameters.config['bedgraph_lower_limit'] == 'auto':
                bottom_limit = 'd'
            else:
                bottom_limit = float(self.parameters.config['bedgraph_lower_limit'][i])

            if self.parameters.config['bedgraph_axis_tics_step'] == 0:
                step = 0
            else:
                step = self.parameters.config['bedgraph_axis_tics_step'][i]

            dictonary = dict(information=self.parameters.config['bedgraph_label'][i],
                             file_name=self.parameters.config['bedgraph'][i], upper_limit=upper_limit,
                             bottom_limit=bottom_limit, step=step)
            adjustable_parameters.append(dictonary)

        global_parameters = dict(start=start, end=end, lvl_of_tics=lvl_of_tics, exons=all_exons, scaffold=scaffold,
                                 pdf=pdf)

        preparated_data = dict(global_parameters=global_parameters, adjustable_parameters=adjustable_parameters,
                               config=self.parameters.config, palette=self.parameters.palette,
                               paths_to_fonts=self.parameters.paths_to_fonts)
        if self.maxmax == 1:
            self.parameters.config['bedgraph_upper_limit'] = ['max']
        return (preparated_data)


class Bedgraph_loader(Loader):
    def __init__(self, preparated_data):
        super().__init__(preparated_data)

    def extract_data(self):
        self.coverages = []
        for i in range(len(self.adjustable_parameters)):
            self.path_to_bedgraph = self.adjustable_parameters[i]['file_name']
            with open(self.path_to_bedgraph) as bedgraph:
                pos_of_start = []
                length_of_interval = []
                coverage_on_interval = []

                for line in bedgraph:
                    line = line.strip()
                    line = line.split('\t')

                    if line[0] == self.global_parameters['scaffold'] and ((
                                                                                  self.global_parameters['end'] >= int(
                                                                              line[1]) >= self.global_parameters[
                                                                                      'start']) or (
                                                                                  self.global_parameters[
                                                                                      'start'] <= int(line[2]) <=
                                                                                  self.global_parameters['end'])):

                        if int(line[1]) < self.global_parameters['start']:
                            line[1] = self.global_parameters['start']
                        if int(line[2]) > self.global_parameters['end']:
                            line[2] = self.global_parameters['end']
                        pos_of_start.append(int(line[1]) - self.global_parameters['start'])
                        length_of_interval.append(int(line[2]) - int(line[1]))
                        coverage_on_interval.append((line[3]))

            coverage = [0 for i in range(int(self.global_parameters['end']) - int(self.global_parameters['start']))]
            for i in range(len(pos_of_start)):
                for q in range(length_of_interval[i]):
                    coverage[pos_of_start[i] + q] = coverage_on_interval[i]

            coverage_str = ''
            count = 0
            for i in coverage:
                count += 1
                coverage_str += str(i)
                if count != len(coverage):
                    coverage_str += ','

            self.coverage = coverage_str.split(',')
            self.coverages.append(self.coverage)

    def translation_of_coordinate(self):

        start = self.global_parameters['start']
        end = self.global_parameters['end']
        self.start = int(start) - int(start)
        self.end = int(end) - int(start)
        self.exons = self.global_parameters['exons']
        for k in range(len(self.exons)):

            exons_row = self.exons[k]['exons'].split(';')
            exons = ''
            if exons_row[0] != '':
                for j in exons_row:
                    i = j.split('-')
                    if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                            self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):
                        exons += str(self.global_parameters['start']) + '-' + str(i[1])
                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                            self.global_parameters['end']):
                        exons += str(i[0]) + '-' + str(i[1])
                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[0]) < int(
                            self.global_parameters['end']) and int(i[1]) >= int(self.global_parameters['end']):
                        exons += str(i[0]) + '-' + str(self.global_parameters['end'])
                        exons += ';'
                    elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                            self.global_parameters['end']):
                        exons += str(self.global_parameters['start']) + '-' + (str(self.global_parameters['end']))
                        exons += ';'

                new_exon = ''
                for i in exons.split(';'):
                    if i != '':
                        new_exon += i
                        if i != exons.split(';')[::-1][1]:
                            new_exon += ';'
                exons = new_exon

            exons = exons.split(';')
            correct_exons = ''

            for j in exons:
                if j != '':
                    q = j.split('-')
                    s = int(q[0])
                    e = int(q[1])
                    left = s - int(start)
                    right = e - int(start)
                    correct_exons += (str(left) + '-' + str(right))
                    if j != exons[::-1][0]:
                        correct_exons += ';'
            exons = correct_exons

            exons = exons.split(';')
            self.exons[k]['exons'] = exons

    def hide_introns(self):
        self.introns = []
        line = []
        for i in range(len(self.coverages[0])):
            line.append(i)

        for l in range(len(self.exons)):
            exons = self.exons[l]['exons'][::-1]
            new_line = []
            for j in exons:
                if j != '':
                    j = j.split('-')
                    start = int(j[0])
                    end = int(j[1])
                    new_line += line[start:end]

            new_line = sorted(new_line)

            introns = []
            count = 0
            if new_line == []:
                introns.append(str(self.start) + '-' + str(self.end))

            for k in range(len(new_line)):
                if count == 0 and new_line[0] != self.start:
                    count += 1
                    introns.append(str(self.start) + '-' + str(new_line[0]))

                if k != len(new_line) - 1:
                    if new_line[k] != (new_line[k + 1] - 1):
                        start = int(new_line[k]) + 1
                        end = new_line[k + 1]
                        introns.append(str(start) + '-' + str(end))
                else:
                    if new_line[len(new_line) - 1] != self.end - 1:
                        start = int(new_line[len(new_line) - 1]) + 1
                        end = self.end
                        introns.append(str(start) + '-' + str(end))
            self.introns.append(dict(introns=introns))

        intersect_introns = []
        for i in range(len(self.introns)):
            if i == 0:
                intersect_introns += set(self.introns[i]['introns'])

            intersect_introns = (set(intersect_introns) & set(self.introns[i]['introns']))

            if i == 0:

                adjust_param = []
                for k in self.introns[i]['introns']:
                    k = k.split('-')
                    start = int(k[0])
                    end = int(k[1])
                    for l in range(end - start):
                        adjust_param.append(start)
                        start += 1
                intersect_introns1 = adjust_param

            if i != 0:
                adjust_param = []

                for p in self.introns[i]['introns']:
                    p = p.split('-')
                    start = int(p[0])
                    end = int(p[1])
                    for q in range(end - start):
                        adjust_param.append(start)
                        start += 1

            intersect_introns1 = (set(intersect_introns1) & set(adjust_param))

        intersect_introns4 = []
        for o in intersect_introns1:
            intersect_introns4.append(o)
        intersect_introns2 = intersect_introns4

        intersect_introns3 = []
        for z in range(len(intersect_introns2) - 1):
            if z == 0:
                start = int(intersect_introns2[0])
            if intersect_introns2[z] != intersect_introns2[z + 1] - 1:
                end = int(intersect_introns2[z]) + 1
                intersect_introns3.append(str(start) + '-' + str(end))
                start = int(intersect_introns2[z + 1])
            if z == len(intersect_introns2) - 1 and end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1]) + 1

                intersect_introns3.append(str(start) + '-' + str(end))
            if z == (len(intersect_introns2) - 2) and self.end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1] + 1)
                intersect_introns3.append(str(start) + '-' + str(end))

        intersect_introns = intersect_introns3

        true_exons = []
        new_line = []
        for i in intersect_introns:
            if i != '':
                i = i.split('-')
                start = int(i[0])
                end = int(i[1])
                new_line += line[start:end]
        new_line = sorted(new_line)
        count = 0
        if new_line == []:
            true_exons.append(str(self.start) + '-' + str(self.end))
        for i in range(len(new_line)):
            if count == 0 and new_line[0] != self.start:
                count += 1
                true_exons.append(str(self.start) + '-' + str(new_line[0]))
            if i != len(new_line) - 1:
                if new_line[i] != (new_line[i + 1] - 1):
                    start = int(new_line[i]) + 1
                    end = new_line[i + 1]
                    true_exons.append(str(start) + '-' + str(end))
            else:
                if new_line[len(new_line) - 1] != self.end - 1:
                    start = int(new_line[len(new_line) - 1]) + 1
                    end = self.end
                    true_exons.append(str(start) + '-' + str(end))

        row_coverages = self.coverages
        coverages = []
        for self.coverage in self.coverages:
            coverage = []
            for i in true_exons:
                SandE = i.split('-')
                start = int(SandE[0])
                end = int(SandE[1])
                coverage += self.coverage[start:end]
            self.coverage = coverage
            coverages.append(self.coverage)
        self.coverages = coverages

    def revcomp_transform(self):

        length = (self.end - self.start)

        for k in range(len(self.exons)):

            new_exon = ''
            exons = self.exons[k]['exons']
            for i in exons:
                if i != '':
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])

                    new_start = (length - end)
                    new_end = (length - start)
                    new_exon += str(new_start) + '-' + str(new_end)
                    if i != exons[::-1][0]:
                        new_exon += ';'
            self.exons[k]['exons'] = new_exon.split(';')

        coverages = []
        for coverage in self.coverages:
            coverage = coverage[::-1]
            coverages.append(coverage)
        self.coverages = coverages

    def create_track(self):
        tracks = []
        count1 = 0
        for k in range(len(self.adjustable_parameters)):

            if count1 != (len(self.config['c_bedgraph_tracks'])):
                color = self.config['c_bedgraph_tracks'][count1]
            else:
                count1 = 0
                color = self.config['c_bedgraph_tracks'][count1]
            count1 += 1

            self.coverage = self.coverages[k]

            if self.global_parameters['lvl_of_tics'] == 'auto':
                if self.adjustable_parameters[k]['upper_limit'] == 0:
                    coverage = []
                    for i in self.coverage:
                        coverage.append(float(i))
                    sorted_coverage = sorted(coverage)
                    max = 0

                else:
                    max = (self.adjustable_parameters[k]['upper_limit'])

                if self.adjustable_parameters[k]['step'] == 0:
                    coverage = []
                    for i in self.coverage:
                        coverage.append(float(i))
                    min = self.adjustable_parameters[k]['bottom_limit']
                    if max != 'd' and min != 'd':
                        max = float(max)
                        min = float(min)
                        if (max - min) >= 3:
                            step = 0
                        else:
                            step = 0
                    else:
                        step = 0
                else:
                    step = self.adjustable_parameters[k]['step']

                axis_tics = []
                min = self.adjustable_parameters[k]['bottom_limit']
                count = min
                if step != 0:
                    while count < max:
                        count += step
                        axis_tics.append(count)
                else:
                    if min != 'd':
                        axis_tics = [min]

            else:
                axis_tics = self.global_parameters['lvl_of_tics']
                min = self.adjustable_parameters[k]['bottom_limit']
                max = self.adjustable_parameters[k]['upper_limit']

            pdf = self.global_parameters['pdf']
            coverage = self.coverage
            information = self.adjustable_parameters[k]['information']

            min = min
            max = max
            try:
                min = float(min)
                max = float(max)

            except:
                pass
            input_data = dict(pdf=self.global_parameters['pdf'], color=color, coverage=self.coverage, upper_limit=max,
                              bottom_limit=min, axis_tics=axis_tics, information=information)

            track = drawing.Bedgraph(input_data, self.config, self.palette, self.paths_to_fonts)

            tracks.append(track)

        return (tracks)


class Paired_bedgraph_data_preparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_dict(self):
        start = self.parameters.gtf['start']
        end = self.parameters.gtf['end']
        step = self.parameters.config['bedgraph_axis_tics_step']
        lvl_of_tics = self.parameters.config['bedgraph_axis_tics']
        self.maxmax = 0
        if self.parameters.config['bedgraph_upper_limit'] != 'd':
            if self.parameters.config['bedgraph_upper_limit'][0] == 'max':
                self.parameters.config['bedgraph_upper_limit'] = 'd'
                self.maxmax = 1
        self.minmin = 0
        if self.parameters.config['bedgraph_lower_limit'] != 'd':
            if self.parameters.config['bedgraph_lower_limit'][0] == 'min':
                self.parameters.config['bedgraph_lower_limit'] = 'd'
                self.minmin = 1
        upper_limit = self.parameters.config['bedgraph_upper_limit']
        bottom_limit = self.parameters.config['bedgraph_lower_limit']
        scaffold = self.parameters.gtf['name_of_scaffold']
        all_exons = []
        for l in range(len(self.parameters.gtf['for_transcript'])):
            exons = self.parameters.gtf['for_transcript'][l]['exons']
            dictonary = dict(exons=exons)
            all_exons.append(dictonary)
        pdf = self.parameters.pdf

        adjustable_parameters = []
        if self.parameters.config['bedgraph_label'] == '':
            self.parameters.config['bedgraph_label'] = []
            for q in self.parameters.config['paired_bedgraph']:
                q = q[0]
                q = os.path.abspath(q)
                q = os.path.basename(q)
                q = os.path.splitext(q)[0]
                self.parameters.config['bedgraph_label'].append(q)

        for i in range(len(self.parameters.config['paired_bedgraph'])):

            if self.parameters.config['bedgraph_upper_limit'] == 'auto' and self.maxmax == 0:
                upper_limit = 'd'
            elif self.maxmax == 0:
                upper_limit = float(self.parameters.config['bedgraph_upper_limit'][i])

            if self.parameters.config['bedgraph_lower_limit'] == 'auto' and self.minmin == 0:
                bottom_limit = 'd'
            elif self.minmin == 0:
                bottom_limit = float(self.parameters.config['bedgraph_lower_limit'][i])

            if self.parameters.config['bedgraph_axis_tics_step'] == 0:
                step = 0
            else:
                step = self.parameters.config['bedgraph_axis_tics_step'][i]

            dictonary = dict(information=self.parameters.config['bedgraph_label'][i],
                             file_name=self.parameters.config['paired_bedgraph'][i], upper_limit=upper_limit,
                             bottom_limit=bottom_limit, step=step)
            adjustable_parameters.append(dictonary)

        global_parameters = dict(start=start, end=end, lvl_of_tics=lvl_of_tics, exons=all_exons, scaffold=scaffold,
                                 pdf=pdf)

        preparated_data = dict(global_parameters=global_parameters, adjustable_parameters=adjustable_parameters,
                               config=self.parameters.config, palette=self.parameters.palette,
                               paths_to_fonts=self.parameters.paths_to_fonts)
        if self.maxmax == 1:
            self.parameters.config['bedgraph_upper_limit'] = ['max']
        if self.minmin == 1:
            self.parameters.config['bedgraph_lower_limit'] = ['min']
        return (preparated_data)


class Paired_bedgraph_loader(Loader):
    def __init__(self, preparated_data):
        super().__init__(preparated_data)

    def extract_data(self):
        self.coverages_plus = []
        self.coverages_minus = []
        for i in range(len(self.adjustable_parameters)):
            self.path_to_bedgraph = self.adjustable_parameters[i]['file_name'][0]
            with open(self.path_to_bedgraph) as bedgraph:
                pos_of_start = []
                length_of_interval = []
                coverage_on_interval = []

                for line in bedgraph:
                    line = line.strip()
                    line = line.split('\t')
                    if line[0] == self.global_parameters['scaffold'] and ((int(self.global_parameters['end']) >= int(
                            line[1]) >= int(self.global_parameters['start'])) or (self.global_parameters[
                                                                                      'start'] <= int(line[2]) <=
                                                                                  self.global_parameters['end'])):

                        if int(line[1]) < self.global_parameters['start']:
                            line[1] = self.global_parameters['start']
                        if int(line[2]) > self.global_parameters['end']:
                            line[2] = self.global_parameters['end']
                        pos_of_start.append(int(line[1]) - self.global_parameters['start'])
                        length_of_interval.append(int(line[2]) - int(line[1]))
                        coverage_on_interval.append((line[3]))

            coverage = [0 for i in range(int(self.global_parameters['end']) - int(self.global_parameters['start']))]
            for i in range(len(pos_of_start)):
                for q in range(length_of_interval[i]):
                    coverage[pos_of_start[i] + q] = coverage_on_interval[i]

            coverage_str = ''
            count = 0
            for i in coverage:
                count += 1
                coverage_str += str(i)
                if count != len(coverage):
                    coverage_str += ','

            self.coverage = coverage_str.split(',')
            self.coverages_plus.append(self.coverage)
        for i in range(len(self.adjustable_parameters)):
            self.path_to_bedgraph = self.adjustable_parameters[i]['file_name'][1]
            with open(self.path_to_bedgraph) as bedgraph:
                pos_of_start = []
                length_of_interval = []
                coverage_on_interval = []

                for line in bedgraph:
                    line = line.strip()
                    line = line.split('\t')
                    if line[0] == self.global_parameters['scaffold'] and ((int(self.global_parameters['end']) >= int(
                            line[1]) >= int(self.global_parameters['start'])) or (self.global_parameters[
                                                                                      'start'] <= int(line[2]) <=
                                                                                  self.global_parameters['end'])):

                        if int(line[1]) < self.global_parameters['start']:
                            line[1] = self.global_parameters['start']
                        if int(line[2]) > self.global_parameters['end']:
                            line[2] = self.global_parameters['end']
                        pos_of_start.append(int(line[1]) - self.global_parameters['start'])
                        length_of_interval.append(int(line[2]) - int(line[1]))
                        coverage_on_interval.append((line[3]))

            coverage = [0 for i in range(int(self.global_parameters['end']) - int(self.global_parameters['start']))]
            for i in range(len(pos_of_start)):
                for q in range(length_of_interval[i]):
                    coverage[pos_of_start[i] + q] = coverage_on_interval[i]

            coverage_str = ''
            count = 0
            for i in coverage:
                count += 1
                coverage_str += str(i)
                if count != len(coverage):
                    coverage_str += ','
            self.coverage = coverage_str.split(',')
            self.coverages_minus.append(self.coverage)

    def translation_of_coordinate(self):

        start = self.global_parameters['start']
        end = self.global_parameters['end']
        self.start = int(start) - int(start)
        self.end = int(end) - int(start)
        self.exons = self.global_parameters['exons']
        for k in range(len(self.exons)):

            exons_row = self.exons[k]['exons'].split(';')
            exons = ''
            if exons_row[0] != '':
                for j in exons_row:
                    i = j.split('-')
                    if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                            self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):
                        exons += str(self.global_parameters['start']) + '-' + str(i[1])
                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                            self.global_parameters['end']):
                        exons += str(i[0]) + '-' + str(i[1])
                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[0]) < int(
                            self.global_parameters['end']) and int(i[1]) >= int(self.global_parameters['end']):
                        exons += str(i[0]) + '-' + str(self.global_parameters['end'])
                        exons += ';'
                    elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                            self.global_parameters['end']):
                        exons += str(self.global_parameters['start']) + '-' + (str(self.global_parameters['end']))
                        exons += ';'

                new_exon = ''
                for i in exons.split(';'):
                    if i != '':
                        new_exon += i
                        if i != exons.split(';')[::-1][1]:
                            new_exon += ';'
                exons = new_exon

            exons = exons.split(';')
            correct_exons = ''

            for j in exons:
                if j != '':
                    q = j.split('-')
                    s = int(q[0])
                    e = int(q[1])
                    left = s - int(start)
                    right = e - int(start)
                    correct_exons += (str(left) + '-' + str(right))
                    if j != exons[::-1][0]:
                        correct_exons += ';'
            exons = correct_exons

            exons = exons.split(';')
            self.exons[k]['exons'] = exons

    def hide_introns(self):
        self.introns = []
        line = []
        for i in range(len(self.coverages_plus[0])):
            line.append(i)

        for l in range(len(self.exons)):
            exons = self.exons[l]['exons'][::-1]
            new_line = []
            for j in exons:
                if j != '':
                    j = j.split('-')
                    start = int(j[0])
                    end = int(j[1])
                    new_line += line[start:end]

            new_line = sorted(new_line)

            introns = []
            count = 0
            if new_line == []:
                introns.append(str(self.start) + '-' + str(self.end))

            for k in range(len(new_line)):
                if count == 0 and new_line[0] != self.start:
                    count += 1
                    introns.append(str(self.start) + '-' + str(new_line[0]))

                if k != len(new_line) - 1:
                    if new_line[k] != (new_line[k + 1] - 1):
                        start = int(new_line[k]) + 1
                        end = new_line[k + 1]
                        introns.append(str(start) + '-' + str(end))
                else:
                    if new_line[len(new_line) - 1] != self.end - 1:
                        start = int(new_line[len(new_line) - 1]) + 1
                        end = self.end
                        introns.append(str(start) + '-' + str(end))
            self.introns.append(dict(introns=introns))

        intersect_introns = []
        for i in range(len(self.introns)):
            if i == 0:
                intersect_introns += set(self.introns[i]['introns'])

            intersect_introns = (set(intersect_introns) & set(self.introns[i]['introns']))

            if i == 0:

                adjust_param = []
                for k in self.introns[i]['introns']:
                    k = k.split('-')
                    start = int(k[0])
                    end = int(k[1])
                    for l in range(end - start):
                        adjust_param.append(start)
                        start += 1
                intersect_introns1 = adjust_param

            if i != 0:
                adjust_param = []

                for p in self.introns[i]['introns']:
                    p = p.split('-')
                    start = int(p[0])
                    end = int(p[1])
                    for q in range(end - start):
                        adjust_param.append(start)
                        start += 1

            intersect_introns1 = (set(intersect_introns1) & set(adjust_param))

        intersect_introns4 = []
        for o in intersect_introns1:
            intersect_introns4.append(o)
        intersect_introns2 = intersect_introns4

        intersect_introns3 = []
        for z in range(len(intersect_introns2) - 1):
            if z == 0:
                start = int(intersect_introns2[0])
            if intersect_introns2[z] != intersect_introns2[z + 1] - 1:
                end = int(intersect_introns2[z]) + 1
                intersect_introns3.append(str(start) + '-' + str(end))
                start = int(intersect_introns2[z + 1])
            if z == len(intersect_introns2) - 1 and end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1]) + 1

                intersect_introns3.append(str(start) + '-' + str(end))
            if z == (len(intersect_introns2) - 2) and self.end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1] + 1)
                intersect_introns3.append(str(start) + '-' + str(end))

        intersect_introns = intersect_introns3

        true_exons = []
        new_line = []
        for i in intersect_introns:
            if i != '':
                i = i.split('-')
                start = int(i[0])
                end = int(i[1])
                new_line += line[start:end]
        new_line = sorted(new_line)
        count = 0
        if new_line == []:
            true_exons.append(str(self.start) + '-' + str(self.end))
        for i in range(len(new_line)):
            if count == 0 and new_line[0] != self.start:
                count += 1
                true_exons.append(str(self.start) + '-' + str(new_line[0]))
            if i != len(new_line) - 1:
                if new_line[i] != (new_line[i + 1] - 1):
                    start = int(new_line[i]) + 1
                    end = new_line[i + 1]
                    true_exons.append(str(start) + '-' + str(end))
            else:
                if new_line[len(new_line) - 1] != self.end - 1:
                    start = int(new_line[len(new_line) - 1]) + 1
                    end = self.end
                    true_exons.append(str(start) + '-' + str(end))

        # row_coverages = self.coverages
        coverages = []
        for self.coverage in self.coverages_plus:
            coverage = []
            for i in true_exons:
                SandE = i.split('-')
                start = int(SandE[0])
                end = int(SandE[1])
                coverage += self.coverage[start:end]
            self.coverage = coverage
            coverages.append(self.coverage)
        self.coverages_plus = coverages

        coverages = []
        for self.coverage in self.coverages_minus:
            coverage = []
            for i in true_exons:
                SandE = i.split('-')
                start = int(SandE[0])
                end = int(SandE[1])
                coverage += self.coverage[start:end]
            self.coverage = coverage
            coverages.append(self.coverage)
        self.coverages_minus = coverages

    def revcomp_transform(self):

        length = (self.end - self.start)

        for k in range(len(self.exons)):

            new_exon = ''
            exons = self.exons[k]['exons']
            for i in exons:
                if i != '':
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])

                    new_start = (length - end)
                    new_end = (length - start)
                    new_exon += str(new_start) + '-' + str(new_end)
                    if i != exons[::-1][0]:
                        new_exon += ';'
            self.exons[k]['exons'] = new_exon.split(';')

        coverages = []
        for coverage in self.coverages_plus:
            coverage = coverage[::-1]
            coverages.append(coverage)
        self.coverages_plus = coverages

        coverages = []
        for coverage in self.coverages_minus:
            coverage = coverage[::-1]
            coverages.append(coverage)
        self.coverages_minus = coverages

        minus = self.coverages_minus
        self.coverages_minus = self.coverages_plus
        self.coverages_plus = minus

    def create_track(self):
        tracks = []
        count1 = 0
        for k in range(len(self.adjustable_parameters)):

            if count1 != (len(self.config['c_bedgraph_tracks'])):
                color = self.config['c_bedgraph_tracks'][count1]
            else:
                count1 = 0
                color = self.config['c_bedgraph_tracks'][count1]
            count1 += 1

            self.coverage_plus = self.coverages_plus[k]
            self.coverage_minus = self.coverages_minus[k]

            if self.global_parameters['lvl_of_tics'] == 'auto':
                if self.adjustable_parameters[k]['upper_limit'] == 0:
                    coverage = []
                    for i in self.coverage:
                        coverage.append(float(i))
                    sorted_coverage = sorted(coverage)
                    max = 0

                else:
                    max = self.adjustable_parameters[k]['upper_limit']

                if self.adjustable_parameters[k]['step'] == 0:
                    coverage = []
                    for i in self.coverage:
                        coverage.append(float(i))
                    min = self.adjustable_parameters[k]['bottom_limit']
                    if max != 'd' and min != 'd':
                        if (max - min) >= 3:
                            step = 0
                        else:
                            step = 0
                    else:
                        step = 0
                else:
                    step = self.adjustable_parameters[k]['step']

                axis_tics = []
                min = self.adjustable_parameters[k]['bottom_limit']
                count = min
                if step != 0:
                    while count < max:
                        count += step
                        axis_tics.append(count)
                else:
                    if min != 'd':
                        axis_tics = [min]

            else:
                axis_tics = self.global_parameters['lvl_of_tics']
                min = self.adjustable_parameters[k]['bottom_limit']
                max = self.adjustable_parameters[k]['upper_limit']

            pdf = self.global_parameters['pdf']
            coverage = self.coverage
            information = self.adjustable_parameters[k]['information']

            min = min
            max = max

            input_data = dict(pdf=self.global_parameters['pdf'], color=color, coverage_plus=self.coverage_plus,
                              coverage_minus=self.coverage_minus, upper_limit=max,
                              bottom_limit=min, axis_tics=axis_tics, information=information)

            track = drawing.Paired_bedgraph(input_data, self.config, self.palette, self.paths_to_fonts)

            tracks.append(track)

        return (tracks)


class Aa_seq_data_preparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)
        self.fasta_file = self.parameters.config['fasta_file']

    def create_dict(self):
        adjustable_parameters = []
        for i in self.parameters.gtf['for_transcript']:
            dictonary = dict(exons=i['exons'])
            adjustable_parameters.append(dictonary)

        start = self.parameters.gtf['start']
        end = self.parameters.gtf['end']
        pdf = self.parameters.pdf
        scaffold = self.parameters.gtf['name_of_scaffold']
        add_codon = self.parameters.config['highlight_codon']
        highlight_reading_frame = self.parameters.config['highlight_reading_frame']
        global_parameters = dict(fasta_file=self.fasta_file, scaffold=scaffold, start=start, add_codon=add_codon,
                                 end=end, pdf=pdf, highlight_reading_frame=highlight_reading_frame)

        preparated_data = dict(global_parameters=global_parameters, adjustable_parameters=adjustable_parameters,
                               config=self.parameters.config, palette=self.parameters.palette,
                               paths_to_fonts=self.parameters.paths_to_fonts)

        return (preparated_data)


class Aa_seq_loader(Loader):
    def __init__(self, preparated_data):
        super().__init__(preparated_data)
        self.path_to_fasta = self.global_parameters['fasta_file']

    def extract_data(self):
        name_of_file = self.path_to_fasta
        name_of_scaffold = self.global_parameters['scaffold']
        start = self.global_parameters['start']
        end = self.global_parameters['end']
        sequnce_obj = data_processing.SequenceExtraction(name_of_file, name_of_scaffold, start, end)
        self.sequence = sequnce_obj.extract().upper()

        region_length = (int(end) - int(start))
        if len(self.sequence) != (region_length):
            fasta = name_of_file
            contig_length = 'N/A'
            for i in (SeqIO.parse(name_of_file, 'fasta')):
                if (str(i.id)) == name_of_scaffold:
                    contig_length = (len(str(i.seq)))
            print("Warning: the window of " + name_of_scaffold + ':' + str(start) + '-' + str(
                end) + ' does not fit the contig ' + name_of_scaffold + ' as given in the supplied fasta file. The detected length of ' + name_of_scaffold + ' is ' + str(
                contig_length), file=sys.stderr)

    def translation_of_coordinate(self):
        for k in range(len(self.adjustable_parameters)):
            self.start = int(self.global_parameters['start'])
            self.end = int(self.global_parameters['end'])
            start = self.start
            end = self.end
            self.start = int(start) - int(start)
            self.end = int(end) - int(start)
            exons_row = self.adjustable_parameters[k]['exons'].split(';')
            exons = ''
            if exons_row[0] != '':
                for j in exons_row:
                    i = j.split('-')

                    if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                            self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):

                        exons += str(self.global_parameters['start']) + '-' + str(i[1])

                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                            self.global_parameters['end']):
                        exons += str(i[0]) + '-' + str(i[1])

                        exons += ';'
                    elif int(i[0]) >= int(self.global_parameters['start']) and int(i[0]) < int(
                            self.global_parameters['end']) and int(i[1]) >= int(self.global_parameters['end']):

                        exons += str(i[0]) + '-' + str(self.global_parameters['end'])

                        exons += ';'
                    elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                            self.global_parameters['end']):
                        exons += str(self.global_parameters['start']) + '-' + (str(self.global_parameters['end']))
                        exons += ';'
                new_exon = ''
                for i in exons.split(';'):
                    if i != '':
                        new_exon += i
                        if i != exons.split(';')[::-1][1]:
                            new_exon += ';'
                exons = new_exon
            else:
                self.exons = ''

            exons = exons.split(';')
            correct_exons = ''

            if len(exons) != 0:
                for j in exons:
                    if j != '':
                        i = j.split('-')
                        s = int(i[0])
                        e = int(i[1])
                        left = s - int(start)
                        right = e - int(start)
                        correct_exons += (str(left) + '-' + str(right))
                        if j != exons[::-1][0]:
                            correct_exons += ';'
                self.exons = correct_exons

            self.adjustable_parameters[k]['exons'] = self.exons.split(';')

    def hide_introns(self):

        line = []
        for i in range(self.end - self.start):
            line.append(i)

        for i in range(len(self.adjustable_parameters)):
            exons = self.adjustable_parameters[i]['exons']
            new_line = []
            for j in exons:
                if j != '':
                    j = j.split('-')
                    start = int(j[0])
                    end = int(j[1])
                    new_line += line[start:end]

            new_line = sorted(new_line)

            introns = []
            count = 0
            if new_line == []:
                introns.append(str(self.start) + '-' + str(self.end))

            for k in range(len(new_line)):
                if count == 0 and new_line[0] != self.start:
                    count += 1
                    introns.append(str(self.start) + '-' + str(new_line[0]))

                if k != len(new_line) - 1:
                    if new_line[k] != (new_line[k + 1] - 1):
                        start = int(new_line[k]) + 1
                        end = new_line[k + 1]
                        introns.append(str(start) + '-' + str(end))
                else:
                    if new_line[len(new_line) - 1] != self.end - 1:
                        start = int(new_line[len(new_line) - 1]) + 1
                        end = self.end
                        introns.append(str(start) + '-' + str(end))
            self.adjustable_parameters[i]['introns'] = introns

        intersect_introns = []
        for i in range(len(self.adjustable_parameters)):
            if i == 0:
                intersect_introns += set(self.adjustable_parameters[i]['introns'])

            intersect_introns = (set(intersect_introns) & set(self.adjustable_parameters[i]['introns']))

            if i == 0:

                adjust_param = []
                for k in (self.adjustable_parameters[i]['introns']):
                    k = k.split('-')
                    start = int(k[0])
                    end = int(k[1])
                    for l in range(end - start):
                        adjust_param.append(start)
                        start += 1
                intersect_introns1 = adjust_param

            if i != 0:
                adjust_param = []

                for p in (self.adjustable_parameters[i]['introns']):
                    p = p.split('-')
                    start = int(p[0])
                    end = int(p[1])
                    for q in range(end - start):
                        adjust_param.append(start)
                        start += 1

            intersect_introns1 = (set(intersect_introns1) & set(adjust_param))

        intersect_introns4 = []
        for o in intersect_introns1:
            intersect_introns4.append(o)
        intersect_introns2 = intersect_introns4

        intersect_introns3 = []
        for z in range(len(intersect_introns2) - 1):
            if z == 0:
                start = int(intersect_introns2[0])
            if intersect_introns2[z] != intersect_introns2[z + 1] - 1:
                end = int(intersect_introns2[z]) + 1
                intersect_introns3.append(str(start) + '-' + str(end))
                start = int(intersect_introns2[z + 1])
            if z == len(intersect_introns2) - 1 and end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1]) + 1

                intersect_introns3.append(str(start) + '-' + str(end))
            if z == (len(intersect_introns2) - 2) and self.end != int(intersect_introns2[z + 1]):
                end = int(intersect_introns2[z + 1] + 1)
                intersect_introns3.append(str(start) + '-' + str(end))

        intersect_introns = intersect_introns3

        true_exon = []
        new_line = []
        for i in intersect_introns:
            if i != '':
                i = i.split('-')
                start = int(i[0])
                end = int(i[1])
                new_line += line[start:end]
        new_line = sorted(new_line)
        count = 0
        if new_line == []:
            true_exon.append(str(self.start) + '-' + str(self.end))
        for i in range(len(new_line)):
            if count == 0 and new_line[0] != self.start:
                count += 1
                true_exon.append(str(self.start) + '-' + str(new_line[0]))
            if i != len(new_line) - 1:
                if new_line[i] != (new_line[i + 1] - 1):
                    start = int(new_line[i]) + 1
                    end = new_line[i + 1]
                    true_exon.append(str(start) + '-' + str(end))
            else:
                if new_line[len(new_line) - 1] != self.end - 1:
                    start = int(new_line[len(new_line) - 1]) + 1
                    end = self.end
                    true_exon.append(str(start) + '-' + str(end))

        row_sequence = self.sequence
        sequence = ''
        for i in true_exon:
            SandE = i.split('-')
            start = int(SandE[0])
            end = int(SandE[1])
            sequence += row_sequence[start:end]
        self.sequence = sequence

    def revcomp_transform(self):

        length = len(self.sequence)

        for k in range(len(self.adjustable_parameters)):

            new_exon = ''
            self.exons = self.adjustable_parameters[k]['exons']
            for i in self.exons:
                if i != '':
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])

                    new_start = (length - end)
                    new_end = (length - start)
                    new_exon += str(new_start) + '-' + str(new_end)
                    if i != self.exons[::-1][0]:
                        new_exon += ';'
            self.adjustable_parameters[k]['exons'] = new_exon.split(';')

        self.reverse = 1
        row_sequence = Seq(self.sequence)
        sequence = row_sequence.reverse_complement()
        self.sequence = str(sequence)

    def create_track(self):
        tracks = []
        highlight_reading_frame = self.global_parameters['highlight_reading_frame']
        sequence = self.sequence
        pdf = self.global_parameters['pdf']
        additional_codon = self.global_parameters['add_codon']
        input_data = dict(pdf=pdf, sequence=sequence, additional_codon=additional_codon,
                          highlight_reading_frame=highlight_reading_frame)
        if len(sequence) < 160 or self.config['aa_track_style'] == 'codons':

            track = drawing.Aminoacid_sequence(input_data, self.config, self.palette, self.paths_to_fonts)

        else:
            track = drawing.Codons(input_data, self.config, self.palette, self.paths_to_fonts)

        tracks.append(track)
        return (tracks)


class Genomic_intervals_data_preparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_dict(self):
        start = self.parameters.gtf['start']
        end = self.parameters.gtf['end']
        pdf = self.parameters.pdf
        scaffold = self.parameters.gtf['name_of_scaffold']
        coordinates = []
        count = start
        for i in range(end - start):
            coordinates.append(count)
            count += 1

        all_exons = []
        for i in range(len(self.parameters.gtf['for_transcript'])):
            exons = self.parameters.gtf['for_transcript'][i]['exons']
            dictonary = dict(exons=exons)
            all_exons.append(dictonary)

        adjustable_parameters = []
        for i in range(len(self.parameters.config['genomic_intervals'])):
            regions = self.parameters.config['genomic_intervals'][i]
            if self.parameters.config['genomic_intervals_label'] != '':
                information = self.parameters.config['genomic_intervals_label'][i]
            else:
                information = ''
            dictonary = dict(coordinates=coordinates, regions=regions, information=information)
            adjustable_parameters.append(dictonary)

        global_parameters = dict(scaffold=scaffold, exons=all_exons, start=start, end=end, pdf=pdf)

        preparated_data = dict(global_parameters=global_parameters, adjustable_parameters=adjustable_parameters,
                               config=self.parameters.config, palette=self.parameters.palette,
                               paths_to_fonts=self.parameters.paths_to_fonts)

        return (preparated_data)


class Genomic_intervals_data_preparation(Data_preparation):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_dict(self):
        start = self.parameters.gtf['start']
        end = self.parameters.gtf['end']
        pdf = self.parameters.pdf
        scaffold = self.parameters.gtf['name_of_scaffold']
        coordinates = []
        count = start
        for i in range(end - start):
            coordinates.append(count)
            count += 1

        all_exons = []
        for i in range(len(self.parameters.gtf['for_transcript'])):
            exons = self.parameters.gtf['for_transcript'][i]['exons']
            dictonary = dict(exons=exons)
            all_exons.append(dictonary)

        adjustable_parameters = []
        for i in range(len(self.parameters.config['genomic_intervals'])):
            regions = self.parameters.config['genomic_intervals'][i]
            if self.parameters.config['genomic_intervals_label'] != '':
                information = self.parameters.config['genomic_intervals_label'][i]
            else:
                information = ''
            dictonary = dict(coordinates=coordinates, regions=regions, information=information)
            adjustable_parameters.append(dictonary)

        global_parameters = dict(scaffold=scaffold, exons=all_exons, start=start, end=end, pdf=pdf)

        preparated_data = dict(global_parameters=global_parameters, adjustable_parameters=adjustable_parameters,
                               config=self.parameters.config, palette=self.parameters.palette,
                               paths_to_fonts=self.parameters.paths_to_fonts)

        return (preparated_data)


class Genomic_intervals_loader(Loader):
    def __init__(self, preparated_data):
        super().__init__(preparated_data)

    def extract_data(self):
        pass

    def translation_of_coordinate(self):
        if self.adjustable_parameters != []:

            self.start = int(self.global_parameters['start'])
            self.end = int(self.global_parameters['end'])
            start = self.start
            end = self.end
            self.start = self.start = self.start
            self.end = self.end - self.start

            self.exons = self.global_parameters['exons']
            for k in range(len(self.exons)):

                exons_row = self.exons[k]['exons'].split(';')
                exons = ''
                if exons_row != '':
                    for j in exons_row:
                        i = j.split('-')
                        if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                                self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):
                            exons += str(self.global_parameters['start']) + '-' + str(i[1])
                            exons += ';'
                        elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                                self.global_parameters['end']):
                            exons += str(i[0]) + '-' + str(i[1])
                            exons += ';'
                        elif int(i[0]) >= int(self.global_parameters['start']) and int(i[0]) < int(
                                self.global_parameters['end']) and int(i[1]) >= int(self.global_parameters['end']):
                            exons += str(i[0]) + '-' + str(self.global_parameters['end'])
                            exons += ';'
                        elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                                self.global_parameters['end']):
                            exons += str(self.global_parameters['start']) + '-' + (str(self.global_parameters['end']))
                            exons += ';'

                    new_exon = ''
                    for i in exons.split(';'):
                        if i != '':
                            new_exon += i
                            if i != exons.split(';')[::-1][1]:
                                new_exon += ';'
                    exons = new_exon

                exons = exons.split(';')
                correct_exons = ''

                for j in exons:
                    if j != '':
                        q = j.split('-')
                        s = int(q[0])
                        e = int(q[1])
                        left = s - int(start)
                        right = e - int(start)
                        correct_exons += (str(left) + '-' + str(right))
                        if j != exons[::-1][0]:
                            correct_exons += ';'
                exons = correct_exons

                exons = exons.split(';')
                self.exons[k]['exons'] = exons

            for k in range(len(self.adjustable_parameters)):

                regions_row = self.adjustable_parameters[k]['regions']
                regions = ''
                if regions_row[0] != '':
                    for j in regions_row:
                        i = j.split('-')

                        if int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) > int(
                                self.global_parameters['start']) and int(i[1]) <= int(self.global_parameters['end']):

                            regions += str(self.global_parameters['start']) + '-' + str(i[1])

                            regions += ';'
                        elif int(i[0]) >= int(self.global_parameters['start']) and int(i[1]) <= int(
                                self.global_parameters['end']):
                            regions += str(i[0]) + '-' + str(i[1])

                            regions += ';'
                        elif int(i[0]) >= int(self.global_parameters['start']) and int(i[0]) < int(
                                self.global_parameters['end']) and int(i[1]) >= int(self.global_parameters['end']):

                            regions += str(i[0]) + '-' + str(self.global_parameters['end'])

                            regions += ';'
                        elif int(i[0]) <= int(self.global_parameters['start']) and int(i[1]) >= int(
                                self.global_parameters['end']):
                            regions += str(self.global_parameters['start']) + '-' + (str(self.global_parameters['end']))
                            regions += ';'
                    new_regions = ''
                    for i in regions.split(';'):
                        if i != '':
                            new_regions += i
                            if i != regions.split(';')[::-1][1]:
                                new_regions += ';'
                    regions = new_regions
                else:
                    self.regions = ''

                regions = regions.split(';')
                correct_regions = ''

                if len(regions) != 0:
                    for j in regions:
                        if j != '':
                            i = j.split('-')
                            s = int(i[0])
                            e = int(i[1])
                            left = s - int(start)
                            right = e - int(start)
                            correct_regions += (str(left) + '-' + str(right))
                            if j != regions[::-1][0]:
                                correct_regions += ';'
                    self.regions = correct_regions

                self.adjustable_parameters[k]['regions'] = self.regions.split(';')

    def hide_introns(self):
        if self.adjustable_parameters != []:

            self.introns = []
            line = []
            for i in range(len(self.coverages[0])):
                line.append(i)

            for l in range(len(self.exons)):
                exons = self.exons[l]['exons'][::-1]
                new_line = []
                for j in exons:
                    if j != '':
                        j = j.split('-')
                        start = int(j[0])
                        end = int(j[1])
                        new_line += line[start:end]

                new_line = sorted(new_line)

                introns = []
                count = 0
                if new_line == []:
                    introns.append(str(self.start) + '-' + str(self.end))

                for k in range(len(new_line)):
                    if count == 0 and new_line[0] != self.start:
                        count += 1
                        introns.append(str(self.start) + '-' + str(new_line[0]))

                    if k != len(new_line) - 1:
                        if new_line[k] != (new_line[k + 1] - 1):
                            start = int(new_line[k]) + 1
                            end = new_line[k + 1]
                            introns.append(str(start) + '-' + str(end))
                    else:
                        if new_line[len(new_line) - 1] != end - 1:
                            start = int(new_line[len(new_line) - 1]) + 1
                            end = self.end
                            introns.append(str(start) + '-' + str(end))
                self.introns.append(dict(introns=introns))

            intersect_introns = []
            for i in range(len(self.introns)):
                if i == 0:
                    intersect_introns += set(self.introns[i]['introns'])

                intersect_introns = (set(intersect_introns) & set(self.introns[i]['introns']))

            true_exons = []
            new_line = []
            for i in intersect_introns:
                if i != '':
                    i = i.split('-')
                    start = int(i[0])
                    end = int(i[1])
                    new_line += line[start:end]
            new_line = sorted(new_line)
            count = 0
            if new_line == []:
                true_exons.append(str(self.start) + '-' + str(self.end))
            for i in range(len(new_line)):
                if count == 0 and new_line[0] != self.start:
                    count += 1
                    true_exons.append(str(self.start) + '-' + str(new_line[0]))
                if i != len(new_line) - 1:
                    if new_line[i] != (new_line[i + 1] - 1):
                        start = int(new_line[i]) + 1
                        end = new_line[i + 1]
                        true_exons.append(str(start) + '-' + str(end))
                else:
                    if new_line[len(new_line) - 1] != self.end - 1:
                        start = int(new_line[len(new_line) - 1]) + 1
                        end = self.end
                        true_exons.append(str(start) + '-' + str(end))

            self.exons = true_exons

            for k in range(len(self.adjustable_parameters)):
                row_coordinates = self.coordinates
                coordinates = []

                for i in self.exons:
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])
                    coordinates += row_coordinates[start:end]
                self.adjustable_parameters[k]['coordinates'] = coordinates

                exist_exon = []

                for i in range(self.end - self.start):
                    exist_exon.append(1)

                for i in self.exons:
                    SandE = i.split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])
                    start_for_end = []
                    count = start

                    for j in range((end - start)):
                        start_for_end.append(count)
                        count += 1

                    for j in start_for_end:
                        exist_exon[j] = 0

                count = 1
                for i in range(len(exist_exon)):
                    if exist_exon[i] == 1:
                        exist_exon[i] = count
                        count += 1

                for i in range(len(exist_exon)):
                    if i != 0:
                        if exist_exon[i] == 0:
                            exist_exon[i] = exist_exon[i - 1]

                exist_exon = exist_exon[::-1]
                exist_exon.append(0)
                exist_exon = exist_exon[::-1]

                regions = self.adjustable_parameters[k]['regions']
                new_regions = ''
                if regions[0] != '':
                    for i in regions:
                        SandE = i.split('-')
                        start = int(SandE[0])
                        end = int(SandE[1])
                        new_start = start - exist_exon[start]
                        new_end = end - exist_exon[end]
                        new_regions += str(new_start) + '-' + str(new_end)
                        if i != regions[::-1][0]:
                            new_regions += ';'

                    self.regions = new_regions
                else:
                    self.regions = ''
                self.adjustable_parameters[k]['regions'] = self.regions.split(';')

    def revcomp_transform(self):
        if self.adjustable_parameters != []:

            length = (self.end - self.start)

            for k in range(len(self.exons)):

                new_exon = ''
                exons = self.exons[k]['exons']
                for i in exons:
                    if i != '':
                        SandE = i.split('-')
                        start = int(SandE[0])
                        end = int(SandE[1])

                        new_start = (length - end)
                        new_end = (length - start)
                        new_exon += str(new_start) + '-' + str(new_end)
                        if i != exons[::-1][0]:
                            new_exon += ';'
                self.exons[k]['exons'] = new_exon.split(';')

            for k in range(len(self.adjustable_parameters)):
                self.reverse = 1
                self.adjustable_parameters[k]['coordinates'] = self.adjustable_parameters[k]['coordinates'][::-1]

                length = len(self.adjustable_parameters[k]['coordinates'])
                new_regions = ''
                self.regions = self.adjustable_parameters[k]['regions']
                for i in self.regions:
                    if i != '':
                        SandE = i.split('-')
                        start = int(SandE[0])
                        end = int(SandE[1])

                        new_start = (length - end)
                        new_end = (length - start)
                        new_regions += str(new_start) + '-' + str(new_end)
                        if i != self.regions[::-1][0]:
                            new_regions += ';'
                self.adjustable_parameters[k]['regions'] = new_regions.split(';')[::-1]

    def create_track(self):
        if self.adjustable_parameters != []:
            tracks = []
            count1 = 0

            for k in range(len(self.adjustable_parameters)):
                self.adjustable_parameters[k]['intersect'] = 1
                if k != (len(self.adjustable_parameters) - 1):
                    regions1 = self.adjustable_parameters[k]['regions'][::-1]
                    regions2 = self.adjustable_parameters[k + 1]['regions'][::-1]
                    if len(regions1) == 1 and len(regions2) == 1:
                        start1 = int(regions1[0].split('-')[0])
                        end1 = int(regions1[0].split('-')[1])
                        start2 = int(regions2[0].split('-')[0])
                        end2 = int(regions2[0].split('-')[1])
                        first = []
                        l = end1
                        for i in range(end1 - start1):
                            first.append(k)
                            l += 1
                        second = []
                        l = end2
                        for i in range(end2 - start2):
                            second.append(k)
                            l += 1
                        result = list(set(first) & set(second))

                        if len(result) < 10:
                            self.adjustable_parameters[k]['intersect'] = 0

            for k in range(len(self.adjustable_parameters)):

                if type(self.config['c_regions']) != list:
                    self.config['c_regions'] = [self.config['c_regions']]

                if count1 != (len(self.config['c_regions'])):
                    color = self.config['c_regions'][count1]
                else:
                    count1 = 0
                    color = self.config['c_regions'][count1]
                count1 += 1

                regions = self.adjustable_parameters[k]['regions'][::-1]
                coordinates = self.adjustable_parameters[k]['coordinates']
                information = self.adjustable_parameters[k]['information']
                intersect = self.adjustable_parameters[k]['intersect']

                pdf = self.global_parameters['pdf']
                input_data = dict(regions=regions, pdf=pdf, information=information, coordinates=coordinates,
                                  color=color, intersect=intersect)
                track = drawing.Genomic_intervals(input_data, self.config, self.palette, self.paths_to_fonts)
                tracks.append(track)
        else:
            tracks = []
        return (tracks)


class Image:
    def __init__(self, vector, parameters, y_0=0):
        self.vector = vector
        self.y_0 = y_0
        self.parameters = parameters

    def go(self):
        try:
            maxs = []
            mins = []
            for i in self.vector:
                y = i.needed_space()
                if type(i) == drawing.Bedgraph:
                    i.draw(0)
                    maxs.append(i.max1)
                if type(i) == drawing.Paired_bedgraph:
                    i.draw(0)
                    maxs.append(i.max1)
                    mins.append(i.min1)
                if type(i) == drawing.Genomic_intervals:
                    if i.intersect == 0:
                        y = -self.parameters.config['gap'] * cm
                if type(i) == drawing.Transcript_structure:
                    if i.new_line == 0:
                        y = -self.parameters.config['gap'] * cm

                self.y_0 += (y + self.parameters.config['gap'] * cm)

            if '.pdf' not in self.parameters.config['output_filename']:
                name_pdf = self.parameters.config['output_filename'] + '.pdf'
            else:
                name_pdf = self.parameters.config['output_filename']

            self.real_pdf = canvas.Canvas(name_pdf, pagesize=(
                self.parameters.config['page_width'] * cm,
                self.y_0 + 2 * self.parameters.config['margin'] * cm))
            # self.real_pdf = canvas.Canvas(self.data['output_filename'], pagesize=(21*cm, 29.7*cm))

            self.y_0 += self.parameters.config['margin'] * cm

            self.exist_info = 0
            self.exist_axis = 0
            count = 0
            for i in self.vector:

                i.input_data['pdf'] = self.real_pdf
                i.canvas = self.real_pdf
                y = i.needed_space()
                if type(i) == drawing.Bedgraph:
                    if self.parameters.config['bedgraph_upper_limit'] != 0:
                        if self.parameters.config['bedgraph_upper_limit'][0] == 'max':
                            maximum = max(maxs)
                            i.config['upper_limit'] = maximum
                            i.max1 = maximum
                if type(i) == drawing.Paired_bedgraph:
                    if self.parameters.config['bedgraph_upper_limit'] != 0:
                        if self.parameters.config['bedgraph_upper_limit'][0] == 'max':
                            maximum = max(maxs)
                            minimum = min(mins)
                            i.config['upper_limit'] = maximum
                            i.config['lower_limit'] = minimum
                            i.max1 = maximum
                            i.min1 = minimum

                i.draw(self.y_0)

                if type(i) == drawing.Nucleotide_sequence:
                    if count == 0:
                        self.y_0_transcript_structure = self.y_0
                    count += 1
                if type(i) == drawing.Title:
                    self.y_0_information = self.y_0
                    self.exist_info = 1
                if type(i) == drawing.Axis_tics:
                    self.y_0_axis_tics = self.y_0
                    self.coordinates = i.coordinates
                    self.exist_axis = 1
                if type(i) == drawing.Transcript_structure:
                    if count == 0:
                        self.y_0_transcript_structure = self.y_0
                    count += 1
                    if i.new_line == 0:
                        y = -self.parameters.config['gap'] * cm

                if type(i) == drawing.Bedgraph:
                    self.coverage = i.coverage_str

                if type(i) == drawing.Genomic_intervals:
                    if i.intersect == 0:
                        y = -self.parameters.config['gap'] * cm

                self.y_0 -= (y + self.parameters.config['gap'] * cm)

        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create the pdf file. There might be a problem with the \'reportlab\' python package.') from error

    def highlighting_the_frame(self):
        try:
            for hf in self.parameters.config['highlight_frame']:
                frame_start = int(hf[0])
                frame_end = int(hf[1])

                coordinates = self.coordinates
                length_of_sequence = len(coordinates)

                width_of_space = ((self.parameters.config['page_width'] * cm - 2 * self.parameters.config[
                    'margin'] * cm) / length_of_sequence) * 0.1
                width_of_box = ((self.parameters.config['page_width'] * cm - 2 * self.parameters.config[
                    'margin'] * cm + width_of_space) - width_of_space * length_of_sequence) / length_of_sequence

                y_upper = self.y_0_transcript_structure
                y_down = self.parameters.config['margin'] * cm + self.parameters.config['gap'] * cm

                while frame_start not in coordinates:
                    if frame_start <= frame_end:
                        frame_start += 1

                start = coordinates.index(frame_start)

                while frame_end not in coordinates:
                    if frame_end >= frame_start:
                        frame_end -= 1
                end = coordinates.index(frame_end)

                x_of_start = self.parameters.config['margin'] * cm + start * (width_of_space + width_of_box)
                x_of_end = self.parameters.config['margin'] * cm + end * (
                        width_of_box + width_of_space) - width_of_space
                color_fill = (methods.hex_to_rgb(self.parameters.palette[self.parameters.config['c_fill_hframe']]))

                aroundcolor = (methods.hex_to_rgb(self.parameters.palette[self.parameters.config['c_stroke_hframe']]))

                self.real_pdf.setStrokeColorRGB(*aroundcolor, self.parameters.config['c_stroke_hframe_alpha'])
                self.real_pdf.setDash(1)
                self.real_pdf.setLineWidth(0.6)
                self.real_pdf.setLineCap(0)

                self.real_pdf.setFillColorRGB(*color_fill, self.parameters.config['c_fill_hframe_alpha'])

                self.real_pdf.rect(x_of_start, y_down, x_of_end - x_of_start, y_upper - y_down, fill=1)

                pdfmetrics.registerFont(TTFont('regular', self.parameters.paths_to_fonts['regular']))
                pdfmetrics.registerFont(TTFont('mono', self.parameters.paths_to_fonts['mono']))

                self.font_size = self.parameters.config['hframe_font_size']
                face = pdfmetrics.getFont('regular').face
                ascent = (face.ascent * self.font_size) / 1000.0
                descent = (face.descent * self.font_size) / 1000.0
                self.height_of_string = (ascent - descent) / 1.
                self.real_pdf.setFont('regular', self.parameters.config['hframe_font_size'])

                color = methods.hex_to_rgb(self.parameters.palette[self.parameters.config['c_text_hframe']])
                self.real_pdf.setStrokeColorRGB(*color, alpha=1)
                self.real_pdf.setFillColorRGB(*color, alpha=1)

                # self.real_pdf.line(x_of_start, y_down, x_of_start, y_down - self.height_of_string/2)
                # self.real_pdf.line(x_of_end, y_down, x_of_end, y_down - self.height_of_string/2)

                y_down -= 1 * self.height_of_string

                # self.real_pdf.drawRightString(x_of_start, y_down, str('{:,}'.format(frame_start)))
                # self.real_pdf.drawString(x_of_end, y_down, str('{:,}'.format(frame_end)))
                self.real_pdf.setFont('regular', self.parameters.config['hframe_font_size'])
                y_down -= 0.1 * self.height_of_string

                x_center = (x_of_end + x_of_start) / 2

                self.real_pdf.drawCentredString(x_center, y_down, hf[2])

        except Exception as error:
            raise methods.Svist4getError(
                'Unable to highlight the requested frame. Please check the entered data (see the program manual).') from error

    def save(self):
        try:
            self.real_pdf.save()
            if '.png' not in self.parameters.config['output_filename']:

                if '.pdf' not in self.parameters.config['output_filename']:
                    name_pdf = self.parameters.config['output_filename'] + '.pdf'
                else:
                    name_pdf = self.parameters.config['output_filename']
                if self.parameters.config['verbose']:
                    print('Image successfully saved to', name_pdf)
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to save the pdf file. There might be a problem with the \'reportlab\' python package.') from error

    def draw(self):

        self.go()
        if self.parameters.config['highlight_frame'] != 0:
            self.highlighting_the_frame()

        self.save()


class Tracks_maker():
    def __init__(self, parameters):
        self.parameters = parameters

    def create_tracks(self):
        pass


class Vgrid_tracks_maker(Tracks_maker):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_tracks(self):
        try:
            obj_with_data = Vgrid_data_preparation(self.parameters)
            preparated_data = obj_with_data.create_dict()
            loader = Vgrid_loader(preparated_data)
            loader.extract_data()
            loader.translation_of_coordinate()

            if self.parameters.config['hide_introns']:
                loader.hide_introns()

            if self.parameters.config['revcomp_transform']:
                loader.revcomp_transform()

            tracks = loader.create_track()

            return (tracks)
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create a Vgrid track. Please check the entered data and the installation of the required modules (see the program manual).') from error


class Title_tracks_maker(Tracks_maker):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_tracks(self):
        try:
            obj_with_data = Title_data_preparation(self.parameters)
            preparated_data = obj_with_data.create_dict()
            loader = Title_loader(preparated_data)
            loader.extract_data()
            loader.translation_of_coordinate()

            if self.parameters.config['hide_introns']:
                loader.hide_introns()

            if self.parameters.config['revcomp_transform']:
                loader.revcomp_transform()

            tracks = loader.create_track()

            return (tracks)
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create a Image label track. Please check the entered data and the installation of the required modules (see the program manual).') from error


class Axis_tics_tracks_maker(Tracks_maker):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_tracks(self):
        try:
            obj_with_data = Axis_tics_data_praparation(self.parameters)
            preparated_data = obj_with_data.create_dict()
            loader = Axis_tics_loader(preparated_data)
            loader.extract_data()
            loader.translation_of_coordinate()

            if self.parameters.config['hide_introns']:
                loader.hide_introns()

            if self.parameters.config['revcomp_transform']:
                loader.revcomp_transform()

            tracks = loader.create_track()

            return (tracks)
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create a Axis tics track. Please check the entered data and the installation of the required modules (see the program manual).') from error


class Nt_seq_tracks_maker(Tracks_maker):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_tracks(self):
        try:
            obj_with_data = Nt_seq_data_preparation(self.parameters)
            preparated_data = obj_with_data.create_dict()
            loader = Nt_seq_loader(preparated_data)
            loader.extract_data()
            loader.translation_of_coordinate()

            if self.parameters.config['hide_introns']:
                loader.hide_introns()

            if self.parameters.config['revcomp_transform']:
                loader.revcomp_transform()

            tracks = loader.create_track()

            return (tracks)
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create a Nucleotide sequence track. Please check the entered data and the installation of the required modules (see the program manual).') from error


class Transcript_struct_tracks_maker(Tracks_maker):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_tracks(self):
        try:
            obj_with_data = Transcript_struct_data_preparation(self.parameters)
            preparated_data = obj_with_data.create_dict()
            loader = Transcript_struct_loader(preparated_data)
            loader.extract_data()
            loader.translation_of_coordinate()

            if self.parameters.config['hide_introns']:
                loader.hide_introns()

            if self.parameters.config['revcomp_transform']:
                loader.revcomp_transform()

            tracks = loader.create_track()
            if self.parameters.config['revcomp_transform']:
                tracks = tracks

            return (tracks)
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create a Transcript structure track. Please check the entered data and the installation of the required modules (see the program manual).') from error


class Bedgraph_tracks_maker(Tracks_maker):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_tracks(self):
        try:
            obj_with_data = Bedgraph_data_preparation(self.parameters)
            preparated_data = obj_with_data.create_dict()
            loader = Bedgraph_loader(preparated_data)
            loader.extract_data()
            loader.translation_of_coordinate()

            if self.parameters.config['hide_introns']:
                loader.hide_introns()

            if self.parameters.config['revcomp_transform']:
                loader.revcomp_transform()

            tracks = loader.create_track()

            return (tracks)
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create a Bedgraph track. Please check the entered data and the installation of the required modules (see the program manual).') from error


class Paired_bedgraph_tracks_maker(Tracks_maker):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_tracks(self):
        try:
            obj_with_data = Paired_bedgraph_data_preparation(self.parameters)
            preparated_data = obj_with_data.create_dict()
            loader = Paired_bedgraph_loader(preparated_data)
            loader.extract_data()
            loader.translation_of_coordinate()

            if self.parameters.config['hide_introns']:
                loader.hide_introns()

            if self.parameters.config['revcomp_transform']:
                loader.revcomp_transform()

            tracks = loader.create_track()

            return (tracks)
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create a Bedgraph track. Please check the entered data and the installation of the required modules (see the program manual).') from error


class Aa_seq_tracks_maker(Tracks_maker):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_tracks(self):
        try:
            obj_with_data = Aa_seq_data_preparation(self.parameters)
            preparated_data = obj_with_data.create_dict()
            loader = Aa_seq_loader(preparated_data)
            loader.extract_data()
            loader.translation_of_coordinate()

            if self.parameters.config['hide_introns']:
                loader.hide_introns()

            if self.parameters.config['revcomp_transform']:
                loader.revcomp_transform()

            tracks = loader.create_track()

            return (tracks)
        except Exception as error:
            raise methods.Svist4getError(
                "Unable to create an Aminoacid sequence track. Please check the entered data and the installation of the required modules (see the program manual). It can be problem with generating or loading fasta sequence file. Please check if 'bedtools' software is installed.") from error


class Genomic_intervals_tracks_maker(Tracks_maker):
    def __init__(self, parameters):
        super().__init__(parameters)

    def create_tracks(self):
        try:
            obj_with_data = Genomic_intervals_data_preparation(self.parameters)
            preparated_data = obj_with_data.create_dict()
            loader = Genomic_intervals_loader(preparated_data)
            loader.extract_data()
            loader.translation_of_coordinate()

            if self.parameters.config['hide_introns']:
                loader.hide_introns()

            if self.parameters.config['revcomp_transform']:
                loader.revcomp_transform()

            tracks = loader.create_track()

            return (tracks)
        except Exception as error:
            raise methods.Svist4getError(
                'Unable to create a Genomic intervals track. Please check the entered data and the installation of the required modules (see the program manual).') from error
