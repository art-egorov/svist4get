from reportlab.pdfgen import canvas
from reportlab.lib.units import cm, mm
from reportlab.lib.colors import black, white
import reportlab.rl_config

reportlab.rl_config.warnOnMissingFontGlyphs = 0
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.lib.pagesizes import letter
from Bio.Seq import Seq
import svist4get.methods as methods
import configs
import math
import numpy as np
import Bio
import statistics


class Track:
    def __init__(self, input_data, config, palette, paths_to_fonts):
        self.config = config
        self.palette = palette
        self.input_data = input_data
        self.paths_to_fonts = paths_to_fonts

    def draw(self):
        pass

    def needed_space(self):
        pass


class Vertical_grid(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        return (0)

    def draw(self, y_upper):
        self.canvas = self.input_data['pdf']
        self.canvas.setStrokeColor(methods.hex_to_rgb(self.palette[self.config['c_vgrid']]),
                                   self.config['c_vgrid_alpha'])
        self.canvas.setDash(0.5, 0.5)
        self.canvas.setLineCap(0)
        self.canvas.setLineWidth(0.45)
        x = self.config['margin'] * cm

        while x <= self.config['page_width'] * cm - self.config['margin'] * cm + 5:
            self.canvas.line(x, y_upper - self.config['gap'] * cm, x,
                             self.config['margin'] * cm + self.config['gap'] * cm)
            x += self.config['vgrid_step'] * cm


class Title(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):

        # pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        # pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))

        self.font_size = self.config['font_size_title']
        face = pdfmetrics.getFont('regular').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_symbol = (ascent - descent) / 1.38

        if self.input_data['additional_info'] != '':

            height = self.height_of_symbol * 2 + 4.5 * self.config['gap'] * cm
        else:
            height = self.height_of_symbol + 2.5 * self.config['gap'] * cm

        return (height)

    def draw(self, y_upper):
        margin = self.config['margin'] * cm

        x = self.config['margin'] * cm

        scaffold = self.input_data['scaffold']
        start = int(self.input_data['start'])
        end = int(self.input_data['end'])
        orientation = self.input_data['strand']

        y = y_upper - self.height_of_symbol

        self.input_data['pdf'].setFillColor(black)
        self.input_data['pdf'].setStrokeColor(black)
        self.input_data['pdf'].setFont('regular', self.font_size)

        first = scaffold + ": " + str('{:,}'.format(start)) + ' - ' + str('{:,}'.format(end)) + ' ' + orientation
        if self.input_data['additional_info'] != '':
            self.input_data['pdf'].drawCentredString(self.config['page_width'] * cm / 2, y,
                                                     self.input_data['additional_info'])
            y = y - self.height_of_symbol - 3 * self.config['gap'] * cm
        self.input_data['pdf'].drawCentredString(self.config['page_width'] * cm / 2, y, first)


class Axis_tics(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        # pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        # pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))
        self.font_size = self.config['gaxis_label_font_size']

        face = pdfmetrics.getFont('regular').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_string_length = (ascent - descent) / 1.38
        self.font_size_length = self.font_size

        self.font_size = self.config['gaxis_tics_font_size']
        face = pdfmetrics.getFont('regular').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_string_tics = (ascent - descent) / 1.38
        self.font_size_tics = self.font_size
        self.height_of_pointer = 0.15 * cm

        if self.config['show_genomic_axis_tics'] != 0:
            self.height = self.height_of_string_length + self.height_of_string_tics + 2 * self.height_of_pointer + 2.5 * \
                          self.config['gap'] * cm
        else:
            self.height = self.height_of_string_length + 1 * self.config['gap'] * cm

        return (self.height)

    def draw(self, y_upper):
        self.step = self.input_data['step']
        self.canvas = self.input_data['pdf']

        height_of_pointer = self.height_of_pointer

        self.coordinates = self.input_data['coordinates']
        coordinates = self.coordinates
        length_of_sequence = len(coordinates)
        y = y_upper - self.height_of_string_length
        self.canvas.setDash(1)
        self.canvas.setStrokeColor(black)
        self.canvas.setLineWidth(0.8)
        self.canvas.setFont('regular', self.font_size_length)
        width_of_string = stringWidth((str(length_of_sequence) + ' bp'), 'regular', self.font_size_length)
        self.canvas.setLineCap(1)
        self.canvas.drawCentredString(self.config['page_width'] * cm / 2, y, ((str(length_of_sequence) + ' bp')))

        y = y_upper - self.height_of_string_length / 2

        self.canvas.line(self.config['margin'] * cm, y,
                         ((self.config['page_width'] / 2 * cm - width_of_string / 2)) - 1.5 * self.config[
                             'gap'] * cm, y)

        self.canvas.line(self.config['margin'] * cm, y,
                         self.config['margin'] * cm + self.height_of_pointer, y + height_of_pointer / 2)
        self.canvas.line(self.config['margin'] * cm, y,
                         self.config['margin'] * cm + self.height_of_pointer, y - height_of_pointer / 2)

        self.canvas.line(self.config['page_width'] * cm - self.config['margin'] * cm, y,
                         ((self.config['page_width'] * cm / 2 + width_of_string / 2)) + 1.5 * self.config[
                             'gap'] * cm, y)
        self.canvas.line(self.config['page_width'] * cm - self.config['margin'] * cm, y,
                         self.config['page_width'] * cm - self.config[
                             'margin'] * cm - self.height_of_pointer, y + height_of_pointer / 2)
        self.canvas.line(self.config['page_width'] * cm - self.config['margin'] * cm, y,
                         self.config['page_width'] * cm - self.config[
                             'margin'] * cm - self.height_of_pointer, y - height_of_pointer / 2)

        if self.config['show_genomic_axis_tics'] != 0:

            y = y_upper - self.height_of_string_length - self.height_of_pointer

            height_of_pointer = self.height_of_pointer
            if self.step == 0:
                self.step = length_of_sequence // 6

            width_of_space = ((self.config['page_width'] * cm - 2 * self.config[
                'margin'] * cm) / length_of_sequence) * 0.1
            width_of_box = ((self.config['page_width'] * cm - 2 * self.config[
                'margin'] * cm + width_of_space) - width_of_space * length_of_sequence) / length_of_sequence

            x = self.config['margin'] * cm
            number_of_cutoffs = (length_of_sequence) // self.step
            length_of_space = self.config['page_width'] * cm - 2 * self.config[
                'margin'] * cm + width_of_space
            length_of_step = (length_of_space) * (self.step / length_of_sequence)

            self.canvas.setFillColor(black)
            self.canvas.setStrokeColor(black)
            self.canvas.setDash(1)
            self.canvas.setLineCap(1)
            self.canvas.setLineWidth(0.8)
            self.canvas.setFont('regular', self.font_size_tics)
            self.canvas.line(x, y, self.config['page_width'] * cm - 1 * self.config['margin'] * cm, y)

            num = int((((length_of_sequence) / 2)))
            x = length_of_space * (num / length_of_sequence) + self.config[
                'margin'] * cm - width_of_space - width_of_box / 2
            deviation = 0
            deviation_num = 0

            while x + deviation < (self.config['page_width'] * cm - self.config['margin'] * cm - 10):
                self.canvas.line(x + deviation, y, x + deviation, y - self.height_of_pointer)
                self.canvas.drawCentredString(x + deviation,
                                              y - self.height_of_pointer - self.height_of_string_tics - self.config[
                                                  'gap'] * cm,
                                              '{:,}'.format(coordinates[num - 1 + deviation_num]))

                deviation += length_of_step
                deviation_num += self.step
            x = length_of_space * (num / length_of_sequence) + self.config[
                'margin'] * cm - width_of_space - width_of_box / 2
            deviation = 0
            deviation_num = 0
            while x - deviation > (self.config['margin'] * cm + 10):
                if deviation != 0:
                    self.canvas.line(x - deviation, y, x - deviation, y - self.height_of_pointer)
                    self.canvas.drawCentredString(x - deviation,
                                                  y - self.height_of_pointer - self.height_of_string_tics - self.config[
                                                      'gap'] * cm,
                                                  '{:,}'.format(coordinates[num - 1 - deviation_num]))
                deviation += length_of_step
                deviation_num += self.step


class Nucleotide_sequence(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        # pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        # pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))

        self.canvas = self.input_data['pdf']
        self.available_width = self.config['page_width'] * cm
        self.sequence = self.input_data['sequence']
        self.length_of_sequence = len(self.sequence)
        self.width_of_space = ((self.available_width - 2 * self.config[
            'margin'] * cm) / self.length_of_sequence) * 0.1
        self.width_of_box = ((self.available_width - 2 * self.config[
            'margin'] * cm + self.width_of_space) - self.width_of_space * self.length_of_sequence) / self.length_of_sequence
        self.width_of_string = stringWidth('A', 'mono', 15)
        self.font_size = 15 * (self.width_of_box / self.width_of_string) * 1.2

        face = pdfmetrics.getFont('mono').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_string = (ascent - descent) / 1.7
        self.deviation = self.height_of_string * 0.15
        self.height_of_box = self.height_of_string * 1.35
        self.canvas.setFont('mono', self.font_size)

        return (self.height_of_box)

    def draw(self, y_upper):

        y = y_upper
        y = y - self.height_of_box
        x = self.config['margin'] * cm
        dx = self.width_of_box + self.width_of_space
        for n in self.sequence:
            if n == "A":
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_A']]),
                                            self.config['c_A_alpha'])
                self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                self.canvas.setFillColor(black)
            if n == "G":
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_G']]),
                                            self.config['c_G_alpha'])
                self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                self.canvas.setFillColor(black)
            if n == "C":
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_C']]),
                                            self.config['c_C_alpha'])
                self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                self.canvas.setFillColor(black)

            if n == "T":
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_T']]),
                                            self.config['c_T_alpha'])
                self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                self.canvas.setFillColor(black)

            elif n != "A" and n != "G" and n != "C" and n != "T":
                self.canvas.setFillColorRGB(0, 0, 0, 0.2)
                self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                self.canvas.setFillColor(black)

            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation, n)
            x = x + dx


class Transcript_structure(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        # pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        # pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))

        self.max_line = self.input_data['max_line']
        self.new_line = self.input_data['new_line']
        self.available_width = self.config['page_width'] * cm
        self.exons = self.input_data['exons']
        self.CDS = self.input_data['CDS']
        self.orientation = self.input_data['orientation']
        self.reverse = self.input_data['reverse']
        self.name = self.input_data['name_of_transcript']

        self.coordinates = self.input_data['coordinates']
        self.delete = self.input_data['delete']
        self.canvas = self.input_data['pdf']

        self.gencode = self.config['triplet_code']

        gencode = ''
        with open(self.gencode) as inp:
            for line in inp:
                gencode += line
        self.gencode = gencode

        self.font_size = float(self.config['transcript_id_font_size'])
        self.width_of_CDS = float(self.config['stroke_width_CDS']) * cm
        self.width_of_exon = float(self.config['stroke_width_exon']) * cm
        face = pdfmetrics.getFont('regular').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_string = (ascent - descent) / 1.38
        self.space = self.config['gap'] * cm

        height = self.height_of_string + self.width_of_CDS + self.space
        return (height)

    def draw(self, y_upper):
        self.canvas.setDash(1)
        margin = self.config['margin'] * cm

        self.length_of_sequence = len(self.coordinates)
        length_of_sequence = self.length_of_sequence

        self.available_width = self.config['page_width'] * cm
        self.width_of_space = ((self.available_width - 2 * self.config[
            'margin'] * cm) / self.length_of_sequence) * 0.1
        width_of_space = self.width_of_space
        self.width_of_box = ((self.available_width - 2 * self.config[
            'margin'] * cm + self.width_of_space) - self.width_of_space * self.length_of_sequence) / self.length_of_sequence
        width_of_box = self.width_of_box

        orientation = self.orientation

        x = margin
        y = y_upper - self.height_of_string

        if self.name == 'No visible annotated transcripts':
            x = self.config['page_width'] * cm / 2
            self.config['transcript_id_font_size'] = 6.5
        else:

            start = self.delete[0].split('-')[1]
            end = self.delete[1].split('-')[0]
            x = (2 * margin + float(start) * (width_of_box + width_of_space) + float(end) * (
                    width_of_box + width_of_space)) / 2
        self.canvas.setFillColor(black)
        self.canvas.setFont('regular', self.config['transcript_id_font_size'])
        width_of_string = stringWidth(self.name, 'regular', self.config['transcript_id_font_size'])

        if x + width_of_string / 2 > self.config['page_width'] * cm - self.config['margin'] * cm:
            self.canvas.drawRightString(self.config['page_width'] * cm - self.config['margin'] * cm, y,
                                        self.name)

        elif x - width_of_string / 2 < self.config['margin'] * cm:
            self.canvas.drawString(self.config['margin'] * cm, y, self.name)
        else:
            self.canvas.drawCentredString(x, y, (self.name))

        y = y - self.space
        y = y - self.width_of_CDS / 2
        self.canvas.setFillColor(black, 1)

        self.transript = self.input_data['transcript']
        x_start = int(self.transript[0]) * (width_of_space + width_of_box) + margin
        x_end = int(self.transript[1]) * (width_of_box + width_of_space) - width_of_space + margin
        self.canvas.setStrokeColor(black, alpha=0.7)
        self.canvas.setLineWidth(1)
        self.canvas.setLineCap(0)
        self.canvas.line(x_start, y, x_end, y)

        self.canvas.setDash(1)
        self.canvas.setLineCap(0)
        self.canvas.setLineWidth(self.width_of_exon)

        self.canvas.setStrokeColor(black, 1)
        if self.exons != '':
            exons_str = self.exons
            for i in exons_str:
                if i != '':
                    i = i.split('-')
                    start = int(i[0])
                    end = int(i[1])

                    x_of_start = margin + (width_of_box + width_of_space) * start
                    x_of_end = margin + (width_of_box + width_of_space) * end - width_of_space
                    while abs(x_of_end - x_of_start) < 0.2:
                        x_of_end += 0.08
                        x_of_start -= 0.08
                    self.canvas.setStrokeColor(white)
                    self.canvas.line(x_of_start, y, x_of_end, y)
                    self.canvas.setStrokeColor(black, 0.7)
                    self.canvas.line(x_of_start, y, x_of_end, y)

        x = margin
        if self.CDS[0] != '':
            self.canvas.setLineWidth(self.width_of_CDS)

            for i in self.CDS:
                i = i.split('-')
                start = int(i[0])
                end = int(i[1])
                x_of_start = margin + (width_of_box + width_of_space) * start
                x_of_end = margin + (width_of_box + width_of_space) * end - width_of_space
                while abs(x_of_end - x_of_start) < 0.2:
                    x_of_end += 0.08
                    x_of_start -= 0.08

                self.canvas.setStrokeColor(white)

                self.canvas.line(x_of_start, y, x_of_end, y)
                self.canvas.setStrokeColor(black, 0.9)

                self.canvas.line(x_of_start, y, x_of_end, y)

        dx = 0.15 * cm
        div = float(self.config['arrow_height']) * cm
        x = margin
        self.canvas.setStrokeColor(black)
        self.canvas.setLineWidth(0.4)
        self.canvas.setDash(1)
        self.canvas.setLineCap(0)
        into = 0
        while x <= self.config['page_width'] * cm:
            if (orientation == '(+)' and self.reverse == 0) or (orientation == '(-)' and self.reverse == 1):

                if x < (margin + (
                        int(self.transript[1]) * (width_of_box + width_of_space)) - width_of_space - 2 * div) and x > (
                        (margin + int(self.transript[0]) * (width_of_space + width_of_box))):

                    if self.exons != '':
                        exons_str = self.exons
                        for i in exons_str:
                            if i != '':
                                i = i.split('-')
                                start = int(i[0])
                                end = int(i[1])

                                x_of_start = margin + (width_of_box + width_of_space) * start
                                x_of_end = margin + (
                                        width_of_box + width_of_space) * end - width_of_space

                                if x >= x_of_start and x <= x_of_end:
                                    into = 1
                    if into == 1:
                        self.canvas.setStrokeColor(white)
                    else:
                        self.canvas.setStrokeColor(black, 0.6)
                    if orientation == '(+)' or orientation == '(-)':
                        self.canvas.setLineCap(1)
                        self.canvas.line(x, y + div, x + div * 2, y)
                        self.canvas.line(x, y - div, x + div * 2, y)
                        self.canvas.setLineCap(0)
                        x += dx
                    into = 0
                else:
                    x += div
                into = 0
            else:

                if x > ((margin + int(self.transript[0]) * (
                        width_of_space + width_of_box)) + 2 * div) and x < (
                        margin + int(self.transript[1]) * (
                        width_of_box + width_of_space) - width_of_space) or x > ((margin + int(self.transript[1]) * (
                        width_of_space + width_of_box)) + 2 * div) and x < (
                        margin + int(self.transript[0]) * (
                        width_of_box + width_of_space) - width_of_space):

                    if self.exons != '':
                        exons_str = self.exons

                        for i in exons_str:
                            if i != '':
                                i = i.split('-')
                                start = int(i[0])
                                end = int(i[1])

                                x_of_start = margin + (width_of_box + width_of_space) * start
                                x_of_end = margin + (
                                        width_of_box + width_of_space) * end - width_of_space

                                if x >= x_of_start and x <= x_of_end:
                                    into = 1
                    if into == 1:
                        self.canvas.setStrokeColor(white)
                    else:
                        self.canvas.setStrokeColor(black, 0.6)
                    if orientation == '(+)' or orientation == '(-)':
                        self.canvas.setLineCap(1)
                        self.canvas.line(x, y + div, x - div * 2, y)
                        self.canvas.line(x, y - div, x - div * 2, y)
                        self.canvas.setLineCap(0)
                    x += dx
                else:
                    x += div
                into = 0

        x = margin
        height = self.config['stroke_width_CDS'] * cm
        height1 = height * 1.4

        for i in range(len(exons_str) - 1):
            r = int((exons_str[i].split('-'))[1])
            l = int((exons_str[i + 1]).split('-')[0])
            r1 = int((exons_str[i].split('-'))[0])
            l1 = int((exons_str[i + 1]).split('-')[1])

            if r == l or r1 == l1:
                if r1 == l1:
                    r = r1
                x = margin + r * (width_of_box + width_of_space) - width_of_space / 2

                if x < (self.available_width - margin) and x > margin:
                    self.canvas.setStrokeColorRGB(1, 1, 1, alpha=1)
                    self.canvas.setLineCap(1)
                    self.canvas.setDash(1)
                    self.canvas.setStrokeColor(white, 1)
                    self.canvas.setLineWidth(0.6)
                    self.canvas.line(x, y + height / 2, x, y - height / 2)

                    self.canvas.setDash(1)
                    self.canvas.setStrokeColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_marks']]), alpha=0.9)
                    self.canvas.setLineWidth(0.4)
                    self.canvas.line(x, y + height1 / 2, x, y - height1 / 2)


'''

class Bedgraph(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        #pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        #pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))

        self.available_width = self.config['page_width'] * cm
        self.coverage_str = self.input_data['coverage']
        self.min1 = self.input_data['bottom_limit']
        self.max1 = self.input_data['upper_limit']
        self.args = self.input_data['axis_tics']
        self.information = self.input_data['information']
        self.canvas = self.input_data['pdf']
        self.color = self.input_data['color']
        height = self.config['bedgraph_track_height'] * cm
        return (height)

    def draw(self, y_upper):
        margin = self.config['margin'] * cm
        gap = self.config['gap'] * cm
        height_for_visualisation_of_coverage = float(self.config['bedgraph_track_height']) * cm
        coverage_str = self.coverage_str

        coverage = []
        for i in coverage_str:
            if self.config['log_scale'] == 0:
                coverage.append(float(i))
            if self.config['log_scale'] == 1:
                coverage.append(math.log2(float(i) + 1))

        num = round((self.available_width - 2 * margin)) / (self.config['bedgraph_column_min_width'] * cm)
        length_of_space = self.available_width - 2 * margin
        length_of_seq = len(coverage)
        width_of_box = length_of_space / num
        current_length = length_of_space / length_of_seq
        coef = length_of_seq / num
        weight = 1
        count = coef
        sum = 0
        k = 0
        new_coverage = []
        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'mean'):
            for i in coverage:
                if (count) // 1 != 0:
                    weight = 1
                    sum += i * weight
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum += i * weight
                    new_coverage.append(sum / coef)
                    sum = 0
                    count = coef
                    sum += i * (1 - weight)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(sum / coef)

            coverage = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0

        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'max'):
            sum = []
            for i in coverage:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(max(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(max(sum))

            coverage = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0


        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'min'):
            sum = []
            for i in coverage:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(min(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(min(sum))

            coverage = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0


        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'median'):
            sum = []
            for i in coverage:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(statistics.median(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(statistics.median(sum))

            coverage = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0

        else:
            length_of_sequence = len(coverage)
            width_of_space = ((self.available_width - 2 * margin) / length_of_sequence) * 0.1
            width_of_box = ((
                                    self.available_width - 2 * margin + width_of_space) - width_of_space * length_of_sequence) / length_of_sequence

        if self.max1 == 'd' or self.min1 == 'd':
            if self.max1 == 'd' and self.min1 != 'd':
                # self.max1 = round(max(coverage_plus) + 0.005, 2)
                self.max1 = abs(self.min1)
            if self.min1 == 'd' and self.max1 != 'd':
                if min(coverage) < 0:
                    self.min1 = self.max1 * -1
                else:
                    self.min1 = 0
            if self.min1 == 'd' and self.max1 == 'd':

                self.max1 = round(max(coverage), self.config['num_of_digits_after_the_dot'])
                self.min1 = round(min(coverage), self.config['num_of_digits_after_the_dot'])
                a = [self.max1, abs(self.min1)]
                self.max1 = max(a)

                if self.min1 < 0:
                    self.min1 = self.max1 * -1
                else:
                    self.min1 = 0

        if self.min1 == 0 and self.max1 == 0:
            if abs(round(max(coverage))) > 0:
                self.max1 = abs(round(max(coverage)))
            else:
                self.max1 = 1


        dx = width_of_box + width_of_space

        y_top = y_upper
        x = margin
        y_bot = y_top - height_for_visualisation_of_coverage
        ratio = height_for_visualisation_of_coverage / (self.max1 - self.min1)
        self.canvas.setLineCap(0)
        count = 0
        y_c = y_bot

        if self.min1 <= 0:
            y_c = y_top - self.max1 * ratio
        elif self.min1 > 0:
            y_c = y_bot


        dx = width_of_box + width_of_space

        self.canvas.setLineCap(0)
        count = 0



        for i in coverage:
            if i >= self.min1 and i <= self.max1:

                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                height = ratio * i
                #if self.min1 > 0:
                 #   height -= ratio * self.min1
                self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)

                self.canvas.setFillColor(black)
                x_c = (x + 0.5 * width_of_box - 0.5)
                # self.canvas.rect(x_c, y_bot, width_of_box * 0.1, 0.5, stroke=0, fill=1)
            elif i < self.min1:
                if self.min1 >= 0:
                    x_c = (x + 0.5 * width_of_box - 0.7)
                    #
                    # self.canvas.setStrokeColorRGB(*methods.hex_to_rgb(self.palette['blue']), 0.9)
                    # self.canvas.rect(x_c, y_bot, width_of_box, width_of_box, stroke=0, fill=1)

                else:
                    self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                    height = ratio * self.min1
                    self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)

                    self.canvas.setFillColorRGB(1, 0, 0, alpha=1)
                    self.canvas.setStrokeColorRGB(1, 0, 0, alpha=1)
                    self.canvas.setLineWidth(0.2)
                    self.canvas.line(x, y_bot - 0.1, x + dx - width_of_space, y_bot - 0.1)
                    x = x + dx
            elif i > self.max1:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                height = ratio * self.max1
                self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)
                x_c = (x + 0.5 * width_of_box - 0.6)
                self.canvas.setFillColorRGB(1, 0, 0, alpha=1)
                self.canvas.setStrokeColorRGB(1, 0, 0, alpha=1)
                self.canvas.setLineWidth(0.2)
                self.canvas.line(x, y_top - 0.1, x + dx - width_of_space, y_top - 0.1)
                # self.canvas.rect(x_c, y_bot, width_of_box, 0.5, stroke=0, fill=1)
            x = x + dx

        points = []
        for i in self.args:
            points.append(float(i))
        if points == []:
            points.append(self.min1)
        for i in points:

            self.canvas.setDash(0.5, 0.5)
            self.canvas.setStrokeColor(black, alpha=0.3)
            self.canvas.setFillColor(black, alpha=self.config['c_bedgraph_label_alpha'])
            self.canvas.setLineWidth(0.05)
            self.canvas.setFont('regular', self.config['bedgraph_tics_font_size'])
            self.font_size = self.config['bedgraph_tics_font_size']
            face = pdfmetrics.getFont('regular').face
            ascent = (face.ascent * self.font_size) / 1000.0
            descent = (face.descent * self.font_size) / 1000.0
            self.height_of_string = (ascent - descent) / 1.38

            self.canvas.setLineWidth(0.3)
            self.canvas.setDash(0.3, 0.1)
            self.canvas.line(margin, y_bot, self.available_width - margin, y_bot)
            self.canvas.setDash(0.3, 0.3)
            self.canvas.line(margin, y_top, self.available_width - margin, y_top)
            self.canvas.drawString(margin, y_bot + self.height_of_string*0.1, str(self.min1))
            self.canvas.drawString(margin, y_top - self.height_of_string*1.1, str(self.max1))

            if i != self.min1:

                point = i - self.min1
                if point <= (self.max1 - self.max1 * 0.08) and point >= self.min1:
                    y_c = y_bot + (y_top - y_bot) * (point / (self.max1 - self.min1))
                    self.canvas.line(margin, y_c, self.available_width - margin, y_c)
                    part_of_range_str = str(point)
                    self.canvas.drawString(margin, y_c - 2.4, part_of_range_str)


            if self.min1 < 0:
                self.canvas.setLineWidth(0.2)
                self.canvas.setStrokeColor(black, 1)
                self.canvas.setDash(1)
                self.canvas.line(margin, y_c, self.available_width - margin, y_c)
            # self.canvas.drawString(margin, y_c, '0')

        self.font_size = self.config['bedgraph_label_font_size']
        face = pdfmetrics.getFont('regular').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_string = (ascent - descent) / 1.38
        self.canvas.setFillColor(black, self.config['c_bedgraph_label_alpha'])
        self.canvas.setFont('regular', self.config['bedgraph_label_font_size'])
        if self.config['bedgraph_label_position'] == 'center':
            self.canvas.drawCentredString(self.available_width / 2, y_top - 1.8 * self.height_of_string, self.information)
        elif self.config['bedgraph_label_position'] == 'left':
            self.canvas.drawString(margin, y_top - 1.8*self.height_of_string, self.information)
        elif self.config['bedgraph_label_position'] == 'right':
            self.canvas.drawRightString(self.available_width - margin, y_top - 1.8*self.height_of_string, self.information)
        else:
            self.canvas.drawCentredString(self.available_width / 2, y_top - 1.8 * self.height_of_string, self.information)

'''

class Bedgraph(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        #pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        #pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))

        self.available_width = self.config['page_width'] * cm
        self.coverage_str = self.input_data['coverage']
        self.min1 = self.input_data['bottom_limit']
        self.max1 = self.input_data['upper_limit']
        self.args = self.input_data['axis_tics']
        self.information = self.input_data['information']
        self.canvas = self.input_data['pdf']
        self.color = self.input_data['color']
        height = self.config['bedgraph_track_height'] * cm

        return (height)

    def draw(self, y_upper):
        margin = self.config['margin'] * cm
        gap = self.config['gap'] * cm
        height_for_visualisation_of_coverage = float(self.config['bedgraph_track_height']) * cm
        coverage_str = self.coverage_str

        coverage = []
        for i in coverage_str:
            if self.config['log_scale'] == 0:
                coverage.append(float(i))
            if self.config['log_scale'] == 1:
                coverage.append(math.log2(float(i) + 1))

        num = round((self.available_width - 2 * margin)) / (self.config['bedgraph_column_min_width'] * cm)
        length_of_space = self.available_width - 2 * margin
        length_of_seq = len(coverage)
        width_of_box = length_of_space / num
        current_length = length_of_space / length_of_seq
        coef = length_of_seq / num
        weight = 1
        count = coef
        sum = 0
        k = 0
        new_coverage = []
        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'mean'):
            for i in coverage:
                if (count) // 1 != 0:
                    weight = 1
                    sum += i * weight
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum += i * weight
                    new_coverage.append(sum / coef)
                    sum = 0
                    count = coef
                    sum += i * (1 - weight)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(sum / coef)

            coverage = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0




        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'max'):
            sum = []
            for i in coverage:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(max(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(max(sum))

            coverage = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0





        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'min'):
            sum = []
            for i in coverage:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(min(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(min(sum))

            coverage = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0


        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'median'):
            sum = []
            for i in coverage:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(statistics.median(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(statistics.median(sum))

            coverage = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0






        else:
            length_of_sequence = len(coverage)
            width_of_space = ((self.available_width - 2 * margin) / length_of_sequence) * 0.1
            width_of_box = ((
                                    self.available_width - 2 * margin + width_of_space) - width_of_space * length_of_sequence) / length_of_sequence

        if self.max1 == 'd' or self.min1 == 'd':
            if self.max1 == 'd' and self.min1 != 'd':
                # self.max1 = round(max(coverage_plus) + 0.005, 2)
                self.max1 = abs(self.min1)
            if self.min1 == 'd' and self.max1 != 'd':
                if min(coverage) < 0:
                    self.min1 = self.max1 * -1
                else:
                    self.min1 = 0
            if self.min1 == 'd' and self.max1 == 'd':

                self.max1 = round(max(coverage), self.config['num_of_digits_after_the_dot'])
                self.min1 = round(min(coverage), self.config['num_of_digits_after_the_dot'])
                if self.max1.is_integer():
                    self.max1 = int(self.max1)
                if self.min1.is_integer():
                    self.min1 = int(self.min1)

                a = [self.max1, abs(self.min1)]
                self.max1 = max(a)

                if self.min1 < 0:
                    self.min1 = self.max1 * -1
                else:
                    self.min1 = 0

        if self.min1 == 0 and self.max1 == 0:
            if abs(round(max(coverage))) > 0:
                self.max1 = abs(round(max(coverage)))
            else:
                self.max1 = 1


        dx = width_of_box + width_of_space

        y_top = y_upper
        x = margin
        y_bot = y_top - height_for_visualisation_of_coverage
        ratio = height_for_visualisation_of_coverage / (self.max1 - self.min1)
        self.canvas.setLineCap(0)
        count = 0
        y_c = y_bot

        if self.min1 <= 0:
            y_c = y_top - self.max1 * ratio
        elif self.min1 > 0:
            y_c = y_bot


        dx = width_of_box + width_of_space

        self.canvas.setLineCap(0)
        count = 0



        for i in coverage:
            if i >= self.min1 and i <= self.max1:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                height = ratio * i
                #if self.min1 > 0:
                 #   height -= ratio * self.min1
                self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)

                self.canvas.setFillColor(black)
                x_c = (x + 0.5 * width_of_box - 0.5)
                # self.canvas.rect(x_c, y_bot, width_of_box * 0.1, 0.5, stroke=0, fill=1)
            elif i < self.min1:
                if self.min1 >= 0:
                    x_c = (x + 0.5 * width_of_box - 0.7)
                    #
                    # self.canvas.setStrokeColorRGB(*methods.hex_to_rgb(self.palette['blue']), 0.9)
                    # self.canvas.rect(x_c, y_bot, width_of_box, width_of_box, stroke=0, fill=1)

                else:
                    self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                    height = ratio * self.min1
                    self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)

                    self.canvas.setFillColorRGB(1, 0, 0, alpha=1)
                    self.canvas.setStrokeColorRGB(1, 0, 0, alpha=1)
                    self.canvas.setLineWidth(0.2)
                    self.canvas.line(x, y_bot - 0.1, x + dx - width_of_space, y_bot - 0.1)
                    x = x + dx
            elif i > self.max1:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                height = ratio * self.max1
                self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)
                x_c = (x + 0.5 * width_of_box - 0.6)
                self.canvas.setFillColorRGB(1, 0, 0, alpha=1)
                self.canvas.setStrokeColorRGB(1, 0, 0, alpha=1)
                self.canvas.setLineWidth(0.2)
                self.canvas.line(x, y_top - 0.1, x + dx - width_of_space, y_top - 0.1)
                # self.canvas.rect(x_c, y_bot, width_of_box, 0.5, stroke=0, fill=1)
            x = x + dx

        points = []
        for i in self.args:
            points.append(float(i))
        if points == []:
            points.append(self.min1)
        for i in points:
            self.canvas.setLineCap(0)
            self.canvas.setDash(0.5, 0.5)
            self.canvas.setStrokeColor(black, alpha=0.3)
            self.canvas.setFillColor(black, alpha=self.config['c_bedgraph_label_alpha'])
            self.canvas.setLineWidth(0.05)
            self.canvas.setFont('regular', self.config['bedgraph_tics_font_size'])
            self.font_size = self.config['bedgraph_tics_font_size']
            face = pdfmetrics.getFont('regular').face
            ascent = (face.ascent * self.font_size) / 1000.0
            descent = (face.descent * self.font_size) / 1000.0
            self.height_of_string = (ascent - descent) / 1.38

            self.canvas.setLineWidth(0.05)
            self.canvas.setStrokeColor(black, alpha=0.3)
            self.canvas.setDash(1)
            self.canvas.line(margin, y_bot - 0.05, self.available_width - margin, y_bot)
            self.canvas.setDash(1)
            self.canvas.setStrokeColor(black, alpha=1)

            self.canvas.setDash(0.5, 0.5)
            self.canvas.setStrokeColor(black, alpha=0.3)
            self.canvas.setLineWidth(0.05)
            self.canvas.line(margin, y_top, self.available_width - margin, y_top)
            self.canvas.drawString(margin, y_top - self.height_of_string * 1.1, str(self.max1))

            if i != self.min1:

                point = i - self.min1
                if point <= (self.max1 - self.max1 * 0.08) and point >= self.min1:
                    y_c = y_bot + (y_top - y_bot) * (point / (self.max1 - self.min1))
                    self.canvas.line(margin, y_c, self.available_width - margin, y_c)
                    part_of_range_str = str(point)
                    self.canvas.drawString(margin, y_c - 2.4, part_of_range_str)


            if self.min1 < 0:
                self.canvas.setLineWidth(0.2)
                self.canvas.setStrokeColor(black, 1)
                self.canvas.setDash(1)
                self.canvas.line(margin, y_c, self.available_width - margin, y_c)
            # self.canvas.drawString(margin, y_c, '0')

        self.font_size = self.config['bedgraph_label_font_size']
        face = pdfmetrics.getFont('regular').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_string = (ascent - descent) / 1.38
        self.canvas.setFillColor(black, self.config['c_bedgraph_label_alpha'])
        self.canvas.setFont('regular', self.config['bedgraph_label_font_size'])
        if self.config['bedgraph_label_position'] == 'center':
            self.canvas.drawCentredString(self.available_width / 2, y_top - 1.8 * self.height_of_string, self.information)
        elif self.config['bedgraph_label_position'] == 'left':
            self.canvas.drawString(margin, y_top - 1.8*self.height_of_string, self.information)
        elif self.config['bedgraph_label_position'] == 'right':
            self.canvas.drawRightString(self.available_width - margin, y_top - 1.8*self.height_of_string, self.information)
        else:
            self.canvas.drawCentredString(self.available_width / 2, y_top - 1.8 * self.height_of_string, self.information)


class Paired_bedgraph(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        #pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        #pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))

        self.available_width = self.config['page_width'] * cm
        self.coverage_str_plus = self.input_data['coverage_plus']
        self.coverage_str_minus = self.input_data['coverage_minus']

        self.min1 = self.input_data['bottom_limit']
        self.max1 = self.input_data['upper_limit']
        self.args = self.input_data['axis_tics']
        self.information = self.input_data['information']
        self.canvas = self.input_data['pdf']
        self.color = self.input_data['color']
        height = self.config['bedgraph_track_height'] * cm

        return (height)

    def draw(self, y_upper):
        margin = self.config['margin'] * cm
        gap = self.config['gap'] * cm
        height_for_visualisation_of_coverage = self.config['bedgraph_track_height'] * cm

        coverage_plus = []
        for i in self.coverage_str_plus:
            if self.config['log_scale'] == 0:
                coverage_plus.append(abs(float(i)))
            else:
                coverage_plus.append(math.log2(abs(float(i)) +1))
        coverage_minus = []
        for i in self.coverage_str_minus:
            if self.config['log_scale'] == 0:
                coverage_minus.append(abs(float(i)))
            else:
                coverage_minus.append(math.log2(abs(float(i))+1))
        num = round((self.available_width - 2 * margin)) / (self.config['bedgraph_column_min_width'] * cm)
        length_of_space = self.available_width - 2 * margin
        length_of_seq = len(coverage_plus)
        width_of_box = length_of_space / num
        current_length = length_of_space / length_of_seq
        coef = length_of_seq / num
        weight = 1
        count = coef
        sum = 0
        k = 0
        new_coverage = []
        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'mean'):
            for i in coverage_plus:
                if (count) // 1 != 0:
                    weight = 1
                    sum += i * weight
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum += i * weight
                    new_coverage.append(sum / coef)
                    sum = 0
                    count = coef
                    sum += i * (1 - weight)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(sum / coef)

            coverage_plus = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0


        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'max'):
            sum = []
            for i in coverage_plus:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(max(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(max(sum))

            coverage_plus = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0



        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'min'):
            sum = []
            for i in coverage_plus:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(min(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(min(sum))

            coverage_plus = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0






        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'median'):
            sum = []
            for i in coverage_plus:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(statistics.median(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(statistics.median(sum))


            coverage_plus = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0
        else:
            length_of_sequence = len(coverage_plus)
            width_of_space = ((self.available_width - 2 * margin) / length_of_sequence) * 0.1
            width_of_box = ((
                                    self.available_width - 2 * margin + width_of_space) - width_of_space * length_of_sequence) / length_of_sequence

        weight = 1
        count = coef
        sum = 0
        k = 0
        new_coverage = []
        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'mean'):
            for i in coverage_minus:
                if (count) // 1 != 0:
                    weight = 1
                    sum += i * weight
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum += i * weight
                    new_coverage.append(sum / coef)
                    sum = 0
                    count = coef
                    sum += i * (1 - weight)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(sum / coef)


            coverage_minus = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0



        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'max'):
            sum = []
            for i in coverage_minus:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(max(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(max(sum))


            coverage_minus = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0



        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'min'):
            sum = []
            for i in coverage_minus:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(min(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(min(sum))


            coverage_minus = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0





        if current_length < self.config['bedgraph_column_min_width'] * cm and (
                self.config['bedgraph_bar'] == 'median'):
            sum = []
            for i in coverage_minus:
                if (count) // 1 != 0:
                    weight = 1
                    sum.append(i)
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum.append(i)
                    new_coverage.append(statistics.median(sum))
                    sum = []
                    count = coef
                    sum.append(i)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(statistics.median(sum))



            coverage_minus = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0
        else:
            length_of_sequence = len(coverage_minus)
            width_of_space = ((self.available_width - 2 * margin) / length_of_sequence) * 0.1
            width_of_box = ((
                                    self.available_width - 2 * margin + width_of_space) - width_of_space * length_of_sequence) / length_of_sequence

        nc = []
        for i in coverage_minus:
            i = i * -1
            nc.append(i)
        coverage_minus = nc

        if self.max1 != 'd' and self.min1 != 'd':
            self.min1 = abs(self.min1)*-1

        if self.max1 == 'd' or self.min1 == 'd':
            if self.max1 == 'd' and self.min1 != 'd':
                # self.max1 = round(max(coverage_plus) + 0.005, 2)
                self.max1 = abs(self.min1)
                self.min1 = abs(self.min1) * -1
            if self.min1 == 'd' and self.max1 != 'd':
                self.min1 = self.max1 * -1

            if self.min1 == 'd' and self.max1 == 'd':
                self.max1 = round(max(coverage_plus), self.config['num_of_digits_after_the_dot'])
                self.min1 = round(min(coverage_minus), self.config['num_of_digits_after_the_dot'])
                if self.max1.is_integer():
                    self.max1 = int(self.max1)
                if self.min1.is_integer():
                    self.min1 = int(self.min1)
                a = [self.max1, abs(self.min1)]
                self.max1 = max(a)
                self.min1 = max(a) * -1

        if self.min1 == 0 and self.max1 == 0:
            self.max1 = 1
            self.min1 = -1

        dx = width_of_box + width_of_space

        y_top = y_upper
        x = margin
        y_bot = y_top - height_for_visualisation_of_coverage

        ratio = height_for_visualisation_of_coverage / (self.max1 - self.min1)
        self.canvas.setLineCap(0)
        count = 0
        y_c = y_bot
        if self.min1 <= 0:
            y_c = y_top - self.max1 * ratio

        # self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette['red']), 0.08)
        # self.canvas.rect(margin, y_c, self.available_width-2*margin, y_top-y_c, stroke = 0, fill =1)
        # self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette['blue']), 0.15)

        # self.canvas.rect(margin, y_bot, self.available_width-2*margin, (y_c-y_bot), stroke = 0, fill =1)

        for i in coverage_plus:
            if i >= self.min1 and i <= self.max1 and i > 0:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                height = ratio * i
                if self.min1 > 0:
                    height -= ratio * self.min1
                self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)
            elif i < self.min1 and i >= 0:
                if self.min1 >= 0:
                    x_c = (x + 0.5 * width_of_box - 0.7)
                    # self.canvas.setStrokeColorRGB(*methods.hex_to_rgb(self.palette['blue']), 0.9)
                    # self.canvas.rect(x_c, y_bot, width_of_box, width_of_box, stroke=0, fill=1)
                else:
                    height = ratio * self.min1
                # self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)
                # self.canvas.setFillColorRGB(1, 0, 0, alpha=1)
                # self.canvas.setStrokeColorRGB(1, 0, 0, alpha=1)
                # self.canvas.setLineWidth(1.2)
                # self.canvas.line(x, y_bot - 0.6, x + dx - width_of_space, y_bot - 0.6)

            elif i > self.max1:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                height = ratio * self.max1
                self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)
                x_c = (x + 0.5 * width_of_box - 0.6)
                self.canvas.setFillColorRGB(1, 0, 0, alpha=1)
                self.canvas.setStrokeColorRGB(1, 0, 0, alpha=1)
                self.canvas.setLineWidth(0.3)
                self.canvas.line(x, y_top - 0.15, x + dx - width_of_space, y_top - 0.15)
                # self.canvas.rect(x_c, y_bot, width_of_box, 0.5, stroke=0, fill=1)
            x = x + dx

        x = margin

        for i in coverage_minus:
            if i >= self.min1 and i <= self.max1 and i <= 0:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                height = ratio * i
                self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)

                # self.canvas.setFillColor(black)
                x_c = (x + 0.5 * width_of_box - 0.5)
                # self.canvas.rect(x_c, y_bot, width_of_box * 0.1, 0.5, stroke=0, fill=1)
            elif i < self.min1 and i < 0 and self.min1 != 0:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                height = ratio * self.min1
                self.canvas.rect(x, y_c, width_of_box, height, stroke=0, fill=1)
                self.canvas.setFillColorRGB(1, 0, 0, alpha=1)
                self.canvas.setStrokeColorRGB(1, 0, 0, alpha=1)
                self.canvas.setLineWidth(0.3)
                self.canvas.line(x, y_bot - 0.15, x + dx - width_of_space, y_bot - 0.15)
            elif i > self.max1:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_bedgraph_alpha'])
                # height = ratio * self.max1
                # self.canvas.rect(x, y_bot, width_of_box, height, stroke=0, fill=1)
                # x_c = (x + 0.5 * width_of_box - 0.6)
                # self.canvas.setFillColorRGB(1, 0, 0, alpha=1)
                # self.canvas.setStrokeColorRGB(1, 0, 0, alpha=1)
                # self.canvas.setLineWidth(1.2)
                # self.canvas.line(x, y_top - 0.6, x + dx - width_of_space, y_top - 0.6)
            x = x + dx

        points = []
        for i in self.args:
            points.append(float(i))

        if points == []:
            points.append(self.min1)
        for i in points:

            self.canvas.setDash(0.5, 0.5)
            self.canvas.setStrokeColor(black, alpha=0.3)
            self.canvas.setFillColor(black, alpha=self.config['c_bedgraph_label_alpha'])
            self.canvas.setLineWidth(0.05)
            self.canvas.setFont('regular', self.config['bedgraph_tics_font_size'])
            self.font_size = self.config['bedgraph_tics_font_size']
            face = pdfmetrics.getFont('regular').face
            ascent = (face.ascent * self.font_size) / 1000.0
            descent = (face.descent * self.font_size) / 1000.0
            self.height_of_string = (ascent - descent) / 1.38

            self.canvas.setLineWidth(0.3)
            self.canvas.setDash(0.3, 0.3)
            self.canvas.line(margin, y_bot, self.available_width - margin, y_bot)
            self.canvas.setDash(0.3, 0.3)
            self.canvas.line(margin, y_top, self.available_width - margin, y_top)
            self.canvas.drawString(margin, y_bot + self.height_of_string * 0.1, str(self.min1))
            self.canvas.drawString(margin, y_top - self.height_of_string * 1.1, str(self.max1))

            # self.canvas.drawRightString(self.available_width-margin, y_bot + 1, str(self.min1))
            # self.canvas.drawRightString(self.available_width - margin, y_top - self.height_of_string, str(self.max1))

            if i != self.min1:

                point = i - self.min1
                if point <= (self.max1 - self.max1 * 0.08) and point >= self.min1:
                    y_c = y_bot + (y_top - y_bot) * (point / (self.max1 - self.min1))
                    self.canvas.line(margin, y_c, self.available_width - margin, y_c)
                    part_of_range_str = str(point)
                    self.canvas.drawString(margin, y_c - 2.4, part_of_range_str)

            if self.min1 < 0:
                self.canvas.setLineWidth(0.2)
                self.canvas.setStrokeColor(black, 1)
                self.canvas.setDash(1)
                self.canvas.line(margin, y_c - 0.1, self.available_width - margin, y_c - 0.1)
            # self.canvas.drawString(margin, y_c, '0')

        self.font_size = self.config['bedgraph_label_font_size']
        face = pdfmetrics.getFont('regular').face
        ascent = (face.ascent * self.font_size*1.1) / 1000.0
        descent = (face.descent * self.font_size*1.1) / 1000.0
        self.height_of_string = (ascent - descent) / 1.38
        self.canvas.setFillColor(black, self.config['c_bedgraph_label_alpha'])
        self.canvas.setFont('regular', self.config['bedgraph_label_font_size'])
        #self.canvas.drawCentredString(self.available_width / 2, y_top - 1.2 * self.height_of_string, self.information)

        if self.config['bedgraph_label_position'] == 'center':
            self.canvas.drawCentredString(self.available_width / 2, y_top - 1.2 * self.height_of_string, self.information)
        elif self.config['bedgraph_label_position'] == 'left':
            self.canvas.drawString(margin, y_top - 1.5*self.height_of_string, self.information)
        elif self.config['bedgraph_label_position'] == 'right':
            self.canvas.drawRightString(self.available_width - margin, y_top - 1.2*self.height_of_string, self.information)
        else:
            self.canvas.drawCentredString(self.available_width / 2, y_top - 1.2 * self.height_of_string, self.information)



'''
class Bedgraph(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        self.available_width = self.config['page_width'] * cm
        self.coverage_str = self.input_data['coverage']
        self.min1 = self.input_data['bottom_limit']
        self.max1 = self.input_data['upper_limit']
        self.args = self.input_data['axis_tics']
        self.information = self.input_data['information']
        self.canvas = self.input_data['pdf']
        self.color = self.input_data['color']
        height = self.config['height_of_bedgraph_track'] * cm

        return (height)

    def draw(self, y_upper):
        margin = self.config['margin'] * cm
        gap = self.config['gap'] * cm
        height_for_visualisation_of_coverage = self.config['height_of_bedgraph_track'] * cm
        coverage_str = self.coverage_str

        coverage = []
        for i in coverage_str:
            coverage.append(float(i))

        num = round((self.available_width - 2 * margin)) / (self.config['minimal_width_of_column'] * cm)
        length_of_space = self.available_width - 2 * margin
        length_of_seq = len(coverage)
        width_of_box = length_of_space / num
        current_length = length_of_space / length_of_seq
        coef = length_of_seq / num
        weight = 1
        count = coef
        sum = 0
        k = 0
        new_coverage = []
        if current_length < self.config['minimal_width_of_column'] * cm and (
                self.config['bedgraph_scale'] == 'auto' or self.config['bedgraph_scale'] == 'small'):
            for i in coverage:
                if (count) // 1 != 0:
                    weight = 1
                    sum += i * weight
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum += i * weight
                    new_coverage.append(sum / coef)
                    sum = 0
                    count = coef
                    sum += i * (1 - weight)
                    count -= (1 - weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(sum / coef)

            coverage = new_coverage
            length_of_sequence = len(new_coverage)
            width_of_box = length_of_space / length_of_sequence
            width_of_space = 0
        else:
            length_of_sequence = len(coverage)
            width_of_space = ((self.available_width - 2 * margin) / length_of_sequence) * 0.1
            width_of_box = ((
                                    self.available_width - 2 * margin + width_of_space) - width_of_space * length_of_sequence) / length_of_sequence

        if self.max1 == 0:
            sorted_coverage = sorted(coverage)
            self.max1 = round(max(coverage) + 0.005, 2)
            if self.max1 == 0:
                self.max1 = 1

        dx = width_of_box + width_of_space

        y_top = y_upper
        x = margin
        y_bot = y_top - height_for_visualisation_of_coverage

        ratio = height_for_visualisation_of_coverage / (self.max1 - self.min1)
        self.canvas.setLineCap(0)
        count = 0
        for i in coverage:
            if i >= self.min1 and i <= self.max1:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_column_alpha'])
                height = ratio * (i - self.min1)
                self.canvas.rect(x, y_bot, width_of_box, height, stroke=0, fill=1)
                self.canvas.setFillColor(black)
                x_c = (x + 0.5 * width_of_box - 0.5)
                self.canvas.rect(x_c, y_bot, width_of_box * 0.1, 0.5, stroke=0, fill=1)
                x = x + dx
            elif i < self.min1:
                x_c = (x + 0.5 * width_of_box - 0.7)
                self.canvas.setStrokeColorRGB(*methods.hex_to_rgb(self.palette['blue']), 0.9)
                self.canvas.rect(x_c, y_bot, width_of_box, width_of_box, stroke=0, fill=1)
                x = x + dx
            elif i > self.max1:
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.color]),
                                            self.config['c_column_alpha'])
                height = ratio * self.max1
                self.canvas.rect(x, y_bot, width_of_box, height, stroke=0, fill=1)
                x_c = (x + 0.5 * width_of_box - 0.6)
                self.canvas.setFillColorRGB(1, 0, 0, alpha=1)
                self.canvas.setStrokeColorRGB(1, 0, 0, alpha=1)
                self.canvas.setLineWidth(1.2)
                self.canvas.line(x, y_top - 0.6, x + dx - width_of_space, y_top - 0.6)
                self.canvas.rect(x_c, y_bot, width_of_box, 0.5, stroke=0, fill=1)
                x = x + dx

        points = []
        for i in self.args:
            points.append(float(i))

        for i in points:

            self.canvas.setDash(0.5, 0.5)
            self.canvas.setStrokeColor(black, alpha=0.3)
            self.canvas.setFillColor(black, alpha=self.config['c_font_size_bedgraph_label_alpha'])
            self.canvas.setLineWidth(0.05)
            self.canvas.setFont('font_medium', self.config['font_size_bedgraph_tics'])
            self.font_size = self.config['font_size_bedgraph_tics']
            face = pdfmetrics.getFont('font_medium').face
            ascent = (face.ascent * self.font_size) / 1000.0
            descent = (face.descent * self.font_size) / 1000.0
            self.height_of_string = (ascent - descent) / 1.6

            if i != self.min1:

                point = i - self.min1
                if point <= (self.max1 - self.max1 * 0.08) and point >= self.min1:
                    y_c = y_bot + (y_top - y_bot) * (point / (self.max1 - self.min1))
                    self.canvas.line(margin, y_c, self.available_width - margin, y_c)
                    part_of_range_str = str(point)
                    self.canvas.drawString(margin, y_c - 2.4, part_of_range_str)

            if i == self.min1:
                self.canvas.setLineWidth(0.3)
                self.canvas.setDash(0.3, 0.1)
                self.canvas.line(margin, y_bot, self.available_width - margin, y_bot)
                self.canvas.setDash(0.3, 0.3)
                self.canvas.line(margin, y_top, self.available_width - margin, y_top)
                self.canvas.drawString(margin, y_bot + 1, str(self.min1))
                self.canvas.drawString(margin, y_top - self.height_of_string, str(self.max1))

        self.font_size = self.config['font_size_bedgraph_label']
        face = pdfmetrics.getFont('font_medium').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_string = (ascent - descent) / 1.7
        self.canvas.setFillColor(black, self.config['c_font_size_bedgraph_label_alpha'])
        self.canvas.setFont('font_regular', self.config['font_size_bedgraph_label'])
        self.canvas.drawCentredString(self.available_width / 2, y_top - 1.8 * self.height_of_string, self.information)
'''

class Aminoacid_sequence(Track):
    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        #pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        #pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))

        self.canvas = self.input_data['pdf']
        self.available_width = self.config['page_width'] * cm
        self.sequence = self.input_data['sequence']
        self.additional_codon = self.input_data['additional_codon']
        self.visability_of_start = self.config['aa_show_start']
        self.gencode = self.config['triplet_code']

        gencode = ''
        with open(self.gencode) as inp:
            for line in inp:
                gencode += line
        self.gencode = gencode

        margin = self.config['margin'] * cm
        gap = self.config['gap'] * cm

        self.length_of_sequence = len(self.sequence)

        self.width_of_space = ((self.available_width - 2 * margin) / self.length_of_sequence) * 0.1
        self.width_of_box = ((
                                     self.available_width - 2 * margin + self.width_of_space - self.length_of_sequence * self.width_of_space) / self.length_of_sequence) * 3 + 2 * self.width_of_space
        self.width_of_box1 = self.width_of_box * 0.8
        self.width_of_string = stringWidth('Cys', 'mono', 15)
        self.font_size = 15 * (self.width_of_box1 / self.width_of_string)

        face = pdfmetrics.getFont('mono').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_string = (ascent - descent) / 1.7
        self.deviation = self.height_of_string * 0.3
        self.height_of_box = self.height_of_string * 1.6
        self.canvas.setFont('mono', self.font_size)
        return (3 * self.height_of_box + 2 * gap)

    def draw(self, y_upper):
        margin = self.config['margin'] * cm
        gap = self.config['gap'] * cm

        self.additional_codon = self.additional_codon.replace('T', 'U')
        # self.gencode = gencode.replace('U', 'T')
        self.gencode = self.gencode.split('\n')
        gendict = {}
        for i in self.gencode:
            if i != '':
                i = i.split('\t')
                val = i[1] + ',' + i[2]
                d = {i[0]: val.split(',')}
                gendict.update(d)

        def get_key(d, value):
            out = []
            for k, v in d.items():
                for i in v:
                    if i == value:
                        out.append(v[0])
                        out.append(k)

            return (out)

        dx = self.width_of_box + self.width_of_space
        y = y_upper - self.height_of_box

        for string in (1, 2, 3):
            small_width_of_box = ((
                                          self.available_width - margin + self.width_of_space) - self.width_of_space * self.length_of_sequence) / self.length_of_sequence
            small_dx = small_width_of_box + self.width_of_space
            x = margin
            if string == 1:
                s = ''
                for i in range(self.length_of_sequence):
                    if len(s) < 3:
                        s += self.sequence[i]
                    if len(s) == 3:
                        s = s.replace('T', 'U')

                        # put stop and start codons in the end of the list with gene code

                        if '!' in get_key(gendict, s) and self.visability_of_start == 1:
                            self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_start_codon']]),
                                                        self.config['c_start_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation, s)


                        elif '#' in get_key(gendict, s):
                            self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_stop_codon']]),
                                                        self.config['c_stop_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation, s)

                        elif s == self.additional_codon:
                            self.canvas.setFillColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_highlight_codon']]),
                                self.config['c_highlight_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation,
                                                          self.additional_codon)


                        else:
                            try:
                                get_key_else = get_key(gendict, s)[1]
                            except:
                                get_key_else = s
                            self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_other_codon']]),
                                                        self.config['c_other_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation,
                                                          get_key_else)

                        x = x + dx
                        s = ''

                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_part_codon']]),
                                            self.config['c_part_codon_alpha'])
                if len(self.sequence) % 3 == 2:
                    self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                    x = x + small_dx
                    self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                if len(self.sequence) % 3 == 1:
                    self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)

            if string == 2:
                y = y - self.height_of_box - gap
                s = ''
                for i in range(self.length_of_sequence):
                    if i == 0:
                        self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_part_codon']]),
                                                    self.config['c_part_codon_alpha'])
                        self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                        x = x + small_dx
                    if len(s) < 3 and i != 0:
                        s += self.sequence[i]
                    if len(s) == 3:
                        s = s.replace('T', 'U')

                        if '!' in get_key(gendict, s) and self.visability_of_start == 1:
                            self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_start_codon']]),
                                                        self.config['c_start_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation, s)


                        elif '#' in get_key(gendict, s):
                            self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_stop_codon']]),
                                                        self.config['c_stop_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation, s)

                        elif s == self.additional_codon:
                            self.canvas.setFillColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_highlight_codon']]),
                                self.config['c_highlight_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation,
                                                          self.additional_codon)


                        else:
                            try:
                                get_key_else = get_key(gendict, s)[1]
                            except:
                                get_key_else = s
                            self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_other_codon']]),
                                                        self.config['c_other_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation,
                                                          get_key_else)

                        x = x + dx
                        s = ''
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_part_codon']]),
                                            self.config['c_part_codon_alpha'])
                if len(self.sequence) % 3 == 0:
                    self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                    x = x + small_dx
                    self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                if len(self.sequence) % 3 == 2:
                    self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                    x = x + small_dx
            if string == 3:
                y = y - self.height_of_box - gap
                s = ''
                for i in range(self.length_of_sequence):
                    if i == 0 or i == 1:
                        self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_part_codon']]),
                                                    self.config['c_part_codon_alpha'])
                        self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                        x = x + small_dx

                    if len(s) < 3 and i != 0 and i != 1:
                        s += self.sequence[i]
                    if len(s) == 3:
                        s = s.replace('T', 'U')

                        if '!' in get_key(gendict, s) and self.visability_of_start == 1:
                            self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_start_codon']]),
                                                        self.config['c_start_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation, s)


                        elif '#' in get_key(gendict, s):
                            self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_stop_codon']]),
                                                        self.config['c_stop_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation, s)

                        elif s == self.additional_codon:
                            self.canvas.setFillColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_highlight_codon']]),
                                self.config['c_highlight_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation,
                                                          self.additional_codon)


                        else:
                            try:
                                get_key_else = get_key(gendict, s)[1]
                            except:
                                get_key_else = s
                            self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_other_codon']]),
                                                        self.config['c_other_codon_alpha'])
                            self.canvas.rect(x, y, self.width_of_box, self.height_of_box, stroke=0, fill=1)
                            self.canvas.setFillColor(black)
                            self.canvas.drawCentredString(x + self.width_of_box / 2, y + self.deviation,
                                                          get_key_else)

                        x = x + dx
                        s = ''
                self.canvas.setFillColorRGB(*methods.hex_to_rgb(self.palette[self.config['c_part_codon']]),
                                            self.config['c_part_codon_alpha'])
                if len(self.sequence) % 3 == 1:
                    self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                    x = x + small_dx
                    self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                if len(self.sequence) % 3 == 0:
                    self.canvas.rect(x, y, small_width_of_box, self.height_of_box, stroke=0, fill=1)
                    x = x + small_dx


class Codons(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        #pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        #pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))

        height = self.config['aa_line_height'] * 3 * 1.25 * cm
        return (height)

    def draw(self, y_upper):
        self.canvas = self.input_data['pdf']
        self.available_width = self.config['page_width'] * cm
        gap = self.config['gap'] * cm
        margin = self.config['margin'] * cm
        self.visability_of_start = self.config['aa_show_start']
        self.additional_codon = self.input_data['additional_codon']
        sequence = self.input_data['sequence']

        length_of_sequence = len(sequence)
        y_top = y_upper
        y_bot = y_upper - self.config['aa_line_height'] * 1.25 * cm

        # 0 for far scale
        self.length_of_sequence = length_of_sequence
        width_of_space = ((self.available_width - 2 * margin) / self.length_of_sequence) * 0.1
        self.width_of_space = width_of_space
        width_of_box = ((
                                self.available_width - 2 * margin + self.width_of_space - self.length_of_sequence * self.width_of_space) / self.length_of_sequence) * 3 + 2 * self.width_of_space

        small_width_of_box = ((
                                      self.available_width - margin + self.width_of_space) - self.width_of_space * self.length_of_sequence) / self.length_of_sequence

        line_height = self.config['aa_line_height'] * cm


        self.gencode = self.config['triplet_code']

        gencode = ''
        with open(self.gencode) as inp:
            for line in inp:
                gencode += line
        self.gencode = gencode

        self.additional_codon = self.additional_codon.replace('T', 'U')

        # self.gencode = gencode.replace('U', 'T')
        self.gencode = gencode.split('\n')
        gendict = {}
        for i in self.gencode:
            if i != '':
                i = i.split('\t')
                val = i[1] + ',' + i[2]
                d = {i[0]: val.split(',')}
                gendict.update(d)

        def get_key(d, value):
            out = []
            for k, v in d.items():
                for i in v:
                    if i == value:
                        out.append(v[0])
                        out.append(k)
            return out

        self.canvas.setDash(1)
        self.canvas.setLineCap(0)

        highlight_frame = int(self.input_data['highlight_reading_frame'])

        for string in (1, 2, 3):

            x = margin
            s = ''
            line_width = self.config['aa_line_width'] * cm
            if string == 1:

                if highlight_frame == 0:
                    x_of_start = self.config['margin'] * cm
                    x_of_end = self.config['page_width'] * cm - self.config['margin'] * cm
                    color_fill = (methods.hex_to_rgb(self.palette[self.config['c_fill_aa_frame']]))

                    aroundcolor = (methods.hex_to_rgb(self.palette[self.config['c_stroke_aa_frame']]))

                    self.canvas.setStrokeColorRGB(*aroundcolor, self.config['c_stroke_aa_frame_alpha'])
                    self.canvas.setDash(1)
                    self.canvas.setLineWidth(0.1)
                    self.canvas.setLineCap(0)

                    self.canvas.setFillColorRGB(*color_fill, self.config['c_fill_aa_frame_alpha'])

                    self.canvas.rect(x_of_start, y_bot, x_of_end - x_of_start, y_top - y_bot, fill=1)

                y_c = (y_top + y_bot) / 2
                self.canvas.setStrokeColorRGB(
                    *methods.hex_to_rgb(self.palette[self.config['c_aa_track_tics_line']]),
                    self.config['c_aa_track_tics_line_alpha'])
                self.canvas.setLineWidth(self.config['aa_line_height'] * cm)
                self.canvas.line(x, y_c, self.available_width - margin, y_c)

                for i in range(length_of_sequence):
                    if len(s) < 3:
                        s += sequence[i]
                    if len(s) == 3:
                        s = s.replace('T', 'U')

                        if '!' in get_key(gendict, s) and self.visability_of_start == 1:
                            self.canvas.setStrokeColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_start_codon']]),
                                self.config['c_start_codon_alpha'])
                            x_c = (x + (x + (width_of_box))) / 2
                            self.canvas.setLineWidth(line_width)
                            self.canvas.line(x_c, y_c + line_height / 2, x_c, y_c - line_height / 2)
                        elif '#' in get_key(gendict, s):
                            self.canvas.setStrokeColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_stop_codon']]),
                                self.config['c_stop_codon_alpha'])
                            x_c = (x + (x + (width_of_box))) / 2
                            self.canvas.setLineWidth(line_width)
                            self.canvas.line(x_c, y_c + line_height / 2, x_c, y_c - line_height / 2)
                        elif s == self.additional_codon:
                            self.canvas.setStrokeColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_highlight_codon']]),
                                self.config['c_highlight_codon_alpha'])
                            x_c = (x + (x + (width_of_box))) / 2
                            self.canvas.setLineWidth(line_width)
                            self.canvas.line(x_c, y_c + line_height / 2, x_c, y_c - line_height / 2)

                        x += width_of_box + width_of_space
                        s = ''

            if string == 2:
                y_top = y_top - self.config['aa_line_height'] * 1.25 * cm
                y_bot = y_bot - self.config['aa_line_height'] * 1.25 * cm

                if highlight_frame == 1:
                    x_of_start = self.config['margin'] * cm
                    x_of_end = self.config['page_width'] * cm - self.config['margin'] * cm
                    color_fill = (methods.hex_to_rgb(self.palette[self.config['c_fill_aa_frame']]))

                    aroundcolor = (methods.hex_to_rgb(self.palette[self.config['c_stroke_aa_frame']]))

                    self.canvas.setStrokeColorRGB(*aroundcolor, self.config['c_stroke_aa_frame_alpha'])
                    self.canvas.setDash(1)
                    self.canvas.setLineWidth(0.1)
                    self.canvas.setLineCap(0)

                    self.canvas.setFillColorRGB(*color_fill, self.config['c_fill_aa_frame_alpha'])

                    self.canvas.rect(x_of_start, y_bot, x_of_end - x_of_start, y_top - y_bot, fill=1)

                y_c = (y_top + y_bot) / 2
                self.canvas.setStrokeColorRGB(
                    *methods.hex_to_rgb(self.palette[self.config['c_aa_track_tics_line']]),
                    self.config['c_aa_track_tics_line_alpha'])
                self.canvas.setLineWidth(self.config['aa_line_height'] * cm)
                self.canvas.line(x, y_c, self.available_width - margin, y_c)

                for i in range(length_of_sequence):
                    if i == 0:
                        x += small_width_of_box + width_of_space

                    if len(s) < 3 and i != 0:
                        s += sequence[i]
                    if len(s) == 3:
                        s = s.replace('T', 'U')

                        if '!' in get_key(gendict, s) and self.visability_of_start == 1:
                            self.canvas.setStrokeColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_start_codon']]),
                                self.config['c_start_codon_alpha'])
                            x_c = (x + (x + (width_of_box))) / 2
                            self.canvas.setLineWidth(line_width)
                            self.canvas.line(x_c, y_c + line_height / 2, x_c, y_c - line_height / 2)
                        elif '#' in get_key(gendict, s):
                            self.canvas.setStrokeColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_stop_codon']]),
                                self.config['c_stop_codon_alpha'])
                            x_c = (x + (x + (width_of_box))) / 2
                            self.canvas.setLineWidth(line_width)
                            self.canvas.line(x_c, y_c + line_height / 2, x_c, y_c - line_height / 2)
                        elif s == self.additional_codon:
                            self.canvas.setStrokeColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_highlight_codon']]),
                                self.config['c_highlight_codon_alpha'])
                            x_c = (x + (x + (width_of_box))) / 2
                            self.canvas.setLineWidth(line_width)
                            self.canvas.line(x_c, y_c + line_height / 2, x_c, y_c - line_height / 2)

                        x += width_of_box + width_of_space

                        s = ''
            if string == 3:
                y_top = y_top - self.config['aa_line_height'] * 1.25 * cm
                y_bot = y_bot - self.config['aa_line_height'] * 1.25 * cm

                if highlight_frame == 2:
                    x_of_start = self.config['margin'] * cm
                    x_of_end = self.config['page_width'] * cm - self.config['margin'] * cm
                    color_fill = (methods.hex_to_rgb(self.palette[self.config['c_fill_aa_frame']]))

                    aroundcolor = (methods.hex_to_rgb(self.palette[self.config['c_stroke_aa_frame']]))

                    self.canvas.setStrokeColorRGB(*aroundcolor, self.config['c_stroke_aa_frame_alpha'])
                    self.canvas.setDash(1)
                    self.canvas.setLineWidth(0.1)
                    self.canvas.setLineCap(0)

                    self.canvas.setFillColorRGB(*color_fill, self.config['c_fill_aa_frame_alpha'])

                    self.canvas.rect(x_of_start, y_bot, x_of_end - x_of_start, y_top - y_bot, fill=1)

                y_c = (y_top + y_bot) / 2
                self.canvas.setStrokeColorRGB(
                    *methods.hex_to_rgb(self.palette[self.config['c_aa_track_tics_line']]),
                    self.config['c_aa_track_tics_line_alpha'])
                self.canvas.setLineWidth(self.config['aa_line_height'] * cm)
                self.canvas.line(x, y_c, self.available_width - margin, y_c)

                for i in range(length_of_sequence):
                    if i == 0 or i == 1:
                        x += small_width_of_box + width_of_space

                    if len(s) < 3 and i != 0 and i != 1:
                        s += sequence[i]
                    if len(s) == 3:
                        s = s.replace('T', 'U')

                        if '!' in get_key(gendict, s) and self.visability_of_start == 1:
                            self.canvas.setStrokeColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_start_codon']]),
                                self.config['c_start_codon_alpha'])
                            x_c = (x + (x + (width_of_box))) / 2
                            self.canvas.setLineWidth(line_width)
                            self.canvas.line(x_c, y_c + line_height / 2, x_c, y_c - line_height / 2)
                        elif '#' in get_key(gendict, s):
                            self.canvas.setStrokeColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_stop_codon']]),
                                self.config['c_stop_codon_alpha'])
                            x_c = (x + (x + (width_of_box))) / 2
                            self.canvas.setLineWidth(line_width)
                            self.canvas.line(x_c, y_c + line_height / 2, x_c, y_c - line_height / 2)
                        elif s == self.additional_codon:
                            self.canvas.setStrokeColorRGB(
                                *methods.hex_to_rgb(self.palette[self.config['c_highlight_codon']]),
                                self.config['c_highlight_codon_alpha'])
                            x_c = (x + (x + (width_of_box))) / 2
                            self.canvas.setLineWidth(line_width)
                            self.canvas.line(x_c, y_c + line_height / 2, x_c, y_c - line_height / 2)

                        x += width_of_box + width_of_space
                        s = ''


class Genomic_intervals(Track):

    def __init__(self, input_data, config, palette, paths_to_fonts):
        super().__init__(input_data, config, palette, paths_to_fonts)

    def needed_space(self):
        # pdfmetrics.registerFont(TTFont('font_medium', self.paths_to_fonts['medium']))
        # pdfmetrics.registerFont(TTFont('font_regular', self.paths_to_fonts['regular']))

        pdfmetrics.registerFont(TTFont('mono', self.paths_to_fonts['mono']))
        pdfmetrics.registerFont(TTFont('regular', self.paths_to_fonts['regular']))

        self.intersect = self.input_data['intersect']
        self.available_width = self.config['page_width'] * cm
        self.regions = self.input_data['regions']
        self.information = self.input_data['information']
        self.canvas = self.input_data['pdf']
        self.coordinates = self.input_data['coordinates']

        self.font_size = float(self.config['font_size_regions_label'])
        self.regions_line_Width = float(self.config['regions_line_width']) * cm
        face = pdfmetrics.getFont('regular').face
        ascent = (face.ascent * self.font_size) / 1000.0
        descent = (face.descent * self.font_size) / 1000.0
        self.height_of_string = (ascent - descent) / 1.38
        self.space = self.config['gap'] * cm
        self.color = self.input_data['color']

        height = self.height_of_string + self.regions_line_Width + self.space / 2

        # height = -2*self.space
        return (height)

    def draw(self, y_upper):
        self.canvas.setDash(1)
        margin = self.config['margin'] * cm

        self.length_of_sequence = len(self.coordinates)
        length_of_sequence = self.length_of_sequence

        self.available_width = self.config['page_width'] * cm
        self.width_of_space = ((self.available_width - 2 * self.config[
            'margin'] * cm) / self.length_of_sequence) * 0.1
        width_of_space = self.width_of_space
        self.width_of_box = ((self.available_width - 2 * self.config[
            'margin'] * cm + self.width_of_space) - self.width_of_space * self.length_of_sequence) / self.length_of_sequence
        width_of_box = self.width_of_box

        if self.regions != '':
            starts = []
            ends = []
            regions_str = self.regions

            for i in regions_str:
                if i != '':
                    i = i.split('-')
                    start = int(i[0])
                    end = int(i[1])

                    x_of_start = margin + (width_of_box + width_of_space) * start
                    x_of_end = margin + (width_of_box + width_of_space) * end - width_of_space
                    starts.append(x_of_start)
                    ends.append(x_of_end)

        x_min = min(starts)
        x_max = max(ends)
        x = (x_min + x_max) / 2

        y = y_upper - self.height_of_string

        self.canvas.setFillColor(black)
        self.canvas.setFont('regular', self.config['font_size_regions_label'])

        self.canvas.drawCentredString(x, y, (self.information))

        y = y - self.regions_line_Width / 2 - self.space / 2
        # y = y - self.regions_line_Width
        self.canvas.setDash(1)
        self.canvas.setLineCap(0)
        self.canvas.setLineWidth(self.regions_line_Width)
        self.canvas.setStrokeColorRGB(*methods.hex_to_rgb(self.palette[self.color]), self.config['c_regions_alpha'])
        if self.regions != '':
            regions_str = self.regions

            for i in regions_str:
                if i != '':
                    i = i.split('-')
                    start = int(i[0])
                    end = int(i[1])
                    x_of_start = margin + (width_of_box + width_of_space) * start
                    x_of_end = margin + (width_of_box + width_of_space) * end - width_of_space

                    self.canvas.line(x_of_start, y, x_of_end, y)

                    # self.canvas.setFont('regular', 13)
                    # self.canvas.drawCentredString(x, y, '▼')

                    self.canvas.setLineWidth(0.63)

                    '''
                    self.canvas.line(x_of_start,y, x_of_start-0.08*cm, y+self.regions_line_Width*0.4)
                    self.canvas.line(x_of_start,y, x_of_start+0.08*cm, y+self.regions_line_Width*0.4)
                    self.canvas.line(x_of_start,y , x_of_start, y+self.regions_line_Width*0.8)

                    
                    self.canvas.line(x_of_start,y, x_of_start-0.1*cm, y+self.regions_line_Width*0.8)
                    self.canvas.line(x_of_start,y, x_of_start+0.1*cm, y+self.regions_line_Width*0.8)
                    self.canvas.line(x_of_start-0.1*cm,y+self.regions_line_Width*0.8, x_of_start+0.1*cm, y+self.regions_line_Width*0.8)
                    

                    
                    self.canvas.line(x_of_start,y, x_of_start, y+self.regions_line_Width*0.8)
                    self.canvas.line(x_of_start,y+self.regions_line_Width*0.8, x_of_start-0.6*cm, y+self.regions_line_Width*0.8)
                    self.canvas.line(x_of_start, y, x_of_start+0.2*self.regions_line_Width, y +  0.3*self.regions_line_Width)
                    self.canvas.line(x_of_start+0.6*cm, y+self.regions_line_Width*0.8, x_of_start+0.6*cm-0.3*self.regions_line_Width, y+self.regions_line_Width*0.8 - 0.2*self.regions_line_Width)
                    '''


'''

class drawing_of_mean_coverage_with_ticks(track):
    def __init__(self, input_data, config = config, palette = palette):
        super().__init__(input_data, config, palette)

    def needed_space(self):
        self.canvas = self.input_data['pdf']
        self.coverage_str = self.input_data['coverage']
        self.available_width = self.config['page_width']*cm
        self.min1 = self.input_data['bottom_limit']
        self.max1 = self.input_data['upper_limit']
        self.information = self.input_data['information']
        self.args = self.input_data['axis_tics']

        height = self.config['height_for_visualisation_of_coverage']*cm
        return(height)

    def draw(self, y_upper):
        margin = self.config['margin']*cm
        height_for_visualisation_of_coverage = self.config['height_for_visualisation_of_coverage']*cm
        gap = self.config['gap']*cm
        coverage = []
        coverage_str = self.coverage_str
        for i in coverage_str:
            coverage.append(float(i))

        num = round((self.available_width - gap))/(0.04*cm)
        length_of_space = self.available_width - 2*margin
        length_of_seq = len(coverage)
        width_of_box = length_of_space / num
        current_length = length_of_space / length_of_seq
        coef = length_of_seq / num
        weight = 1
        count = coef
        sum = 0
        k = 0
        new_coverage = []
        if current_length < 0.04*cm:
            for i in coverage:
                if (count)//1 != 0:
                    weight = 1
                    sum += i*weight
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum += i*weight
                    new_coverage.append(sum/coef)
                    sum = 0
                    count = coef
                    sum += i*(1-weight)
                    count -= (1-weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(sum/coef)
        else:
            new_coverage = coverage
            width_of_box = length_of_space/length_of_seq

        length_of_sequence = len(new_coverage)
        width_of_box = length_of_space/length_of_sequence
        if self.max1 == 0:
            sorted_coverage = sorted(new_coverage)
            self.max1 = round(max(sorted_coverage),2)
            if self.max1 == 0:
                self.max1 = 1

        x = margin
        y_top = y_upper
        y_bot = y_upper - height_for_visualisation_of_coverage
        ratio = height_for_visualisation_of_coverage / (self.max1 - self.min1)
        dx = width_of_box
        self.canvas.setDash(1)

        for i in new_coverage:
            if i >= self.min1 and i <= self.max1:
                self.canvas.setFillColorRGB(*advis4seq_methods.hex_to_rgb(palette[config['c_column']]), self.config['c_column_alpha'])
                height = ratio* (i - self.min1)
                self.canvas.rect(x, y_bot, width_of_box, height, stroke = 0, fill = 1)
                x = x + dx

            elif i < self.min1:
                x = x + dx

            elif i > self.max1:
                self.canvas.setFillColorRGB(*advis4seq_methods.hex_to_rgb(palette[config['c_column']]), self.config['c_column_alpha'])
                height = ratio * self.max1
                self.canvas.rect(x,y_bot, width_of_box, height, stroke = 0, fill = 1)
                self.canvas.setFillColorRGB(1,0,0, alpha=1)
                self.canvas.rect(x, y_top - (height/120) , width_of_box ,(height/120), stroke = 0, fill = 1)
                x = x + dx

        points = []
        for i in self.args:
            points.append(float(i))


        for i in points:

                self.canvas.setDash(0.5,0.5)
                self.canvas.setStrokeColor(black, alpha = 0.3)
                self.canvas.setFillColor(black)
                self.canvas.setLineWidth(0.05)
                self.canvas.setFont('font_medium', 7.5)

                if i != self.min1:

                    point = i - self.min1
                    if point <=(self.max1 - self.max1*0.1)  and point >= (self.min1+ (self.max1*0.05)):
                        y_c = y_bot + (y_top-y_bot)*(point/(self.max1-self.min1))
                        self.canvas.line(margin, y_c, self.available_width - margin, y_c)
                        part_of_range_str = str(point)

                        self.canvas.drawRightString(margin, y_c-2.4, part_of_range_str)

                if  i == self.min1:




                    self.canvas.setStrokeColor(black)
                    self.canvas.setLineWidth(0.3)
                    self.canvas.setDash(0.3, 0.1)
                    self.canvas.line(margin,y_bot, self.available_width-margin , y_bot)
                    self.canvas.setDash(0.3,0.3)
                    self.canvas.line(margin, y_top, self.available_width - margin, y_top)
                    self.canvas.drawRightString(margin, y_bot, str(self.min1))
                    self.canvas.drawRightString(margin, y_top-4.5, str(self.max1))
        self.canvas.setFillColor(black)
        self.canvas.setFont('font_regular', 8)
        self.canvas.drawCentredString(self.available_width/2, (y_top+y_bot)/2 + (y_top-y_bot)/3, self.information)

'''

'''
class drawing_of_mean_coverage_with_step(track):
    def __init__(self, canvas, available_width, coverage_str, information, step = 0, min1=0, max1=0 ):
        self.canvas = canvas
        self.coverage_str = coverage_str
        self.available_width = available_width
        self.min1 = min1
        self.max1 = max1
        self.step = step
        self.information = information
    def needed_space(self):
        height = height_for_visualisation_of_coverege
        return(height)

    def draw(self, y_upper):
        coverage_str = self.coverage_str.split(';')

        coverage = []
        for i in coverage_str:
            coverage.append(float(i))



        num = round((self.available_width - self.config['gap']*cm))/(0.04*cm)
        length_of_space = self.available_width - 2*margin
        length_of_seq = len(coverage)
        width_of_box = length_of_space / num
        current_length = length_of_space / length_of_seq
        coef = length_of_seq / num
        weight = 1
        count = coef
        sum = 0
        k = 0
        new_coverage = []
        if current_length < 0.04*cm:
            for i in coverage:
                if (count)//1 != 0:
                    weight = 1
                    sum += i*weight
                    count -= weight
                    k += 1
                else:
                    weight = count
                    sum += i*weight
                    new_coverage.append(sum/coef)
                    sum = 0
                    count = coef
                    sum += i*(1-weight)
                    count -= (1-weight)
                    k += 1
                if k == length_of_seq:
                    new_coverage.append(sum/coef)
        else:
            new_coverage = coverage
            width_of_box = length_of_space/length_of_seq

        length_of_sequence = len(new_coverage)

        if self.max1 == 0:
            sorted_coverage = sorted(new_coverage)
            self.max1 = sorted_coverage[int(length_of_sequence/1.05)]
            if self.max1 == 0:
                self.max1 = 1
            self.max1 =  round(self.max1, 2)


        if self.step == 0:
            if (self.max1 - self.min1) >= 3:
                self.step = (self.max1 - self.min1) //3
            else:
                self.step =  (self.max1 - self.min1) / 3
                self.step = round(self.step, 2 )

        x = margin
        y_top = y_upper
        y_bot = y_upper - height_for_visualisation_of_coverege
        ratio = height_for_visualisation_of_coverege / (self.max1 - self.min1)
        dx = width_of_box
        self.canvas.setDash(1)

        for i in new_coverage:
            if i >= self.min1 and i <= self.max1:
                self.canvas.setFillColorRGB(*deep_orange_column)
                height = ratio* (i - self.min1)
                self.canvas.rect(x, y_bot, width_of_box, height, stroke = 0, fill = 1)
                x = x + dx

            elif i < self.min1:
                x = x + dx

            elif i > self.max1:
                self.canvas.setFillColorRGB(*deep_orange_column)
                height = ratio * self.max1
                self.canvas.rect(x,y_bot, width_of_box, height, stroke = 0, fill = 1)
                self.canvas.setFillColorRGB(1,0,0, alpha=1)
                self.canvas.rect(x, y_bot , width_of_box , (height/25), stroke = 0, fill = 1)
                x = x + dx



        for i in range(int(((self.max1 - self.min1)/self.step) +1 )):

                self.canvas.setDash(0.5,0.5)
                self.canvas.setStrokeColor(black, alpha = 0.3)
                self.canvas.setFillColor(black)
                self.canvas.setLineWidth(0.05)
                self.canvas.setFont('font_medium', 6)

                if i != 0:

                    part_of_range = self.step*i
                    if part_of_range < (self.max1-self.min1)*0.8:
                        y_c = y_bot + (y_top-y_bot)*(part_of_range/(self.max1-self.min1))
                        self.canvas.line(margin, y_c, self.available_width - margin, y_c)
                        part_of_range_str = str(part_of_range + self.min1)
                        self.canvas.drawRightString(margin, y_c-2.4, part_of_range_str)

                if  i == 0:

                    self.canvas.setStrokeColor(black)
                    self.canvas.setLineWidth(0.3)
                    self.canvas.setDash(0.3, 0.1)
                    self.canvas.line( margin,y_bot, self.available_width-margin , y_bot)
                    self.canvas.setDash(0.3,0.3)
                    self.canvas.line(margin, y_top, self.available_width - margin, y_top)
                    self.canvas.drawRightString(margin, y_bot, str(self.min1))
                    self.canvas.drawRightString(margin, y_top-4, str(self.max1))
        self.canvas.setFillColor(black)
        self.canvas.setFont('fontr_regular', 8)
        self.canvas.drawCentredString(self.available_width/2, (y_top+y_bot)/2 + (y_top-y_bot)/3, self.information)



class drawing_of_coverage_with_step(track):
    def __init__(self, canvas, available_width, coverage_str,  information, step = 0,  min = 0, max = 0):
        self.canvas = canvas
        self.available_width = available_width
        self. step = step
        self.coverage_str = coverage_str
        self.min = min
        self.max = max
        self.information = information

    def needed_space(self):
        height = height_for_visualisation_of_coverege

        return(height)


    def draw(self, y_upper):
        coverage_str = self.coverage_str.split(';')
        self.canvas.setDash(1)
        coverage = []

        for i in coverage_str:
            coverage.append(float(i))

        length_of_sequence = len(coverage)
        width_of_space = ((self.available_width - 2*margin)/length_of_sequence)*0.17
        width_of_box =((self.available_width-2*margin + width_of_space) - width_of_space*length_of_sequence)/length_of_sequence
        dx = width_of_box + width_of_space




        y_top = y_upper
        y_bot = y_top - height_for_visualisation_of_coverege

        self.canvas.setFillColor(black)
        self.canvas.setLineWidth(0.05)
        self.canvas.setFont('fontr_regular', 10)
        self.canvas.saveState()
        self.canvas.translate(self.available_width-margin, (y_bot+y_top)/2)
        self.canvas.rotate(-90)
        self.canvas.drawCentredString(0, 0, self.information)
        self.canvas.restoreState()
  


        x = margin


        if self.max == 0:
            sorted_coverage = sorted(coverage)
            self.max = sorted_coverage[int(length_of_sequence/1.5)]
            if self.max == 0:
                self.max = 1

        if self.step == 0:
            if (self.max - self.min) >= 3:
                self.step = (self.max - self.min) //3
            else:
                self.step =  (self.max - self.min) / 3

        ratio = height_for_visualisation_of_coverege / self.max
        self.canvas.setLineCap(0)


        for i in coverage:
            if i >= self.min and i <= self.max:
                self.canvas.setFillColorRGB(*orange_column)
                height  = ratio * (i-self.min)
                self.canvas.rect(x, y_bot, width_of_box, height, stroke = 0, fill = 1)
                self.canvas.setFillColor(black)
                x_c = (x + 0.5*width_of_box-0.5)
                self.canvas.rect(x_c, y_bot, 1, 1, stroke = 0, fill = 1)
                x = x+dx
            elif i < self.min:
                x_c = (x + 0.5*self.width_of_box - 0.7)
                self.canvas.setStrokeColorRGB(*blue_column)
                self.canvas.rect(x_c, y_bot, 1.4, 1.4, stroke = 0, fill =1)
                x = x + dx
            elif i > self.max:
                self.canvas.setFillColorRGB(*orange_column)
                height = ratio * self.max
                self.canvas.rect(x, y_bot,width_of_box, height, stroke = 0, fill = 1)
                x_c = (x + 0.5*width_of_box - 0.6)
                self.canvas.setFillColorRGB(1,0,0, alpha = 1)
                self.canvas.setStrokeColorRGB(1,0,0, alpha = 1)
                self.canvas.setLineWidth(1.2)
                self.canvas.line(x, y_top-0.6, x + dx - width_of_space, y_top-0.6)
                self.canvas.rect(x_c, y_bot, 1.4,1.4, stroke = 0 , fill =1 )
                x = x+dx

        for i in range(int(((self.max - self.min)//self.step) +1 )):

                self.canvas.setDash(0.5,0.5)
                self.canvas.setStrokeColor(black, alpha = 0.3)
                self.canvas.setFillColor(black)
                self.canvas.setLineWidth(0.05)
                self.canvas.setFont('font_medium', 6)

                if i != 0:

                    part_of_range = self.step*i
                    if part_of_range < (self.max-self.min)*0.8:
                        y_c = y_bot + (y_top-y_bot)*(part_of_range/(self.max-self.min))
                        self.canvas.line(margin, y_c, self.available_width - margin, y_c)
                        part_of_range_str = str(part_of_range + self.min)
                        self.canvas.drawRightString(margin, y_c-2.4, part_of_range_str)

                if  i == 0:

                    self.canvas.setStrokeColor(black)
                    self.canvas.setLineWidth(0.3)
                    self.canvas.setDash(0.3, 0.1)
                    self.canvas.line( margin,y_bot, self.available_width-margin , y_bot)
                    self.canvas.setDash(0.3,0.3)
                    self.canvas.line(margin, y_top, self.available_width - margin, y_top)
                    self.canvas.drawRightString(margin, y_bot, str(self.min))
                    self.canvas.drawRightString(margin, y_top-4, str(self.max))
        self.canvas.setFillColor(black)
        self.canvas.setFont('fontr_regular', 10)
        self.canvas.drawCentredString(self.available_width/2, (y_top+y_bot)/2 + (y_top-y_bot)/3, self.information)

'''
'''
class drawing_of_multi_coverage(track):
    def __init__(self, canvas, available_width, intervals, step, max, min, reverse, *coverages):
        self.canvas = canvas
        self.available_width = available_width
        self.intervals = intervals
        self.step = step
        self.max = max
        self.min = min
        self.row_coverages = coverages
        self.reverse = reverse


    def needed_space(self):
        height = height_for_visualisation_of_coverege
        return(height)

    def draw(self, y_upper):

        length_of_sequence = 0
        if self.intervals != '':
            intervals_str = self.intervals.split(';')
            for i in intervals_str:
                SandE = str(i).split('-')
                start = int(SandE[0])
                end = int(SandE[1])
                length_of_sequence += (end - start)
        else:
            length_of_sequence = len(self.row_coverages[0])



        coverages = [[0]*length_of_sequence for i in range(len(self.row_coverages))]

        for i in range(len(self.row_coverages)):

            row_coverage_str = self.row_coverages[i].split(',')
            coverage_str = []
            if self.intervals != '':
                for j in intervals_str:
                    SandE = str(j).split('-')
                    start = int(SandE[0])
                    end = int(SandE[1])
                    coverage_str += row_coverage_str[start:end]
            else:
                coverage_str = row_coverage_str

            coverage = []
            for k in coverage_str:
                coverage.append(float(k))

            coverages[i] = coverage
            if self.reverse == 1:
                coverages[i] = coverages[i][::-1]


        y_top = y_upper
        y_bot = y_top - height_for_visualisation_of_coverege
        x = margin
        width_of_space = ((self.available_width - 2*margin)/length_of_sequence)*0.17
        width_of_box =((self.available_width-2*margin + width_of_space) - width_of_space*length_of_sequence)/length_of_sequence
        dx = width_of_box + width_of_space
        y = y_bot

        ratio = height_for_visualisation_of_coverege /(self.max - self.min)
        self.canvas.setLineCap(0)

        self.canvas.setLineWidth(0.4)
        self.canvas.setDash(1)
        for j in range(len(self.row_coverages)):
            x = margin
            count = 0
            self.canvas.setStrokeColorRGB(*deep_orange_column)
            for i in coverages[j]:


                if  j ==0:
                    self.canvas.setStrokeColorRGB(0.48, 0, 0.79, 1)
                else:
                    self.canvas.setStrokeColorRGB(1, 0.65, 0, 1)

                try:
                    next = coverages[j][count+1]
                except:
                    next = i
                count +=1
                if i >= self.min and i <= self.max:
                    height1  = ratio * (i-self.min)
                    self.canvas.line(x-width_of_space/2, y+height1, x+width_of_box+width_of_space/2, y+height1)
                    if next <= self.max and next >= self.min:
                        height2 = ratio*(next-self.min)
                    elif next <= self.min:
                        height2 = 0
                    elif next > self.min:
                        height2 = height_for_visualisation_of_coverege

                    self.canvas.line(x+width_of_space/2 + width_of_box, y+height1, x+width_of_box + width_of_space/2,y+ height2)

                    x = x+dx

                elif i < self.min:
                    height1 = 0

                    x_c = (x + 0.5*width_of_box - 0.7)
                    self.canvas.setFillColorRGB(*blue_column)
                    self.canvas.rect(x_c, y_bot, 1.4, 1.4, stroke = 0, fill =1)

                    if next <= self.max and next >= self.min:
                        height2 = ratio*(next-self.min)
                    elif next <= self.min:
                        height2 = 0
                    elif next > self.min:
                        height2 = height_for_visualisation_of_coverege

                    self.canvas.line(x+width_of_space/2 + width_of_box, y+height1, x+width_of_box + width_of_space/2,y+ height2)

                    x = x+dx



                elif i > self.max:
                    height1 = height_for_visualisation_of_coverege


                    #self.canvas.line(x-width_of_space/2, y+height1, x+width_of_box+width_of_space/2, y+height1)

                    x_c = (x + 0.5*width_of_box - 0.7)
                    self.canvas.setFillColorRGB(1,0,0,1)
                    #self.canvas.rect(x_c, y_top-1.4, 1.4, 1.4, stroke = 0, fill =1)


                    if next <= self.max and next >= self.min:
                        height2 = ratio*(next-self.min)
                    elif next <= self.min:
                        height2 = 0
                    elif next > self.min:
                        height2 = height_for_visualisation_of_coverege
                    self.canvas.line(x+width_of_space/2 + width_of_box, y+height1, x+width_of_box + width_of_space/2,y+ height2)

                    x = x+dx

        for i in range(int(((self.max - self.min)//self.step) +1 )):

                self.canvas.setDash(0.5,0.5)
                self.canvas.setStrokeColor(black, alpha = 0.3)
                self.canvas.setFillColor(black)
                self.canvas.setLineWidth(0.05)
                self.canvas.setFont('font_medium', 6)

                if i != 0:

                    part_of_range = self.step*i
                    if part_of_range < (self.max-self.min)*0.8:
                        y_c = y_bot + (y_top-y_bot)*(part_of_range/(self.max-self.min))
                        self.canvas.line(margin, y_c, self.available_width - margin, y_c)
                        part_of_range_str = str(part_of_range + self.min)
                        self.canvas.drawRightString(margin, y_c-2.4, part_of_range_str)

                if  i == 0:

                    self.canvas.setStrokeColor(black)
                    self.canvas.setLineWidth(0.3)
                    self.canvas.setDash(0.3, 0.1)
                    self.canvas.line( margin,y_bot, self.available_width-margin , y_bot)
                    self.canvas.setDash(0.3,0.3)
                    self.canvas.line(margin, y_top, self.available_width - margin, y_top)
                    self.canvas.drawRightString(margin, y_bot, str(self.min))
                    self.canvas.drawRightString(margin, y_top-4, str(self.max))

'''

'''

class drawing_of_CDS(track):

    def __init__(self, canvas, config, palette, variable_data):
        super().__init__(canvas, config, palette, variable_data)
    
    def __init__(self, canvas, available_width, sequence, CDS, reverse):
        self.canvas  = canvas
        self.available_width = available_width
        self.sequence = sequence
        self.CDS = CDS
        self.reverse = reverse
        
    def needed_space(self):
        return(0.3*cm)

    def draw(self, y_upper):
        sequence = self.sequence
        sequence = self.sequence


        length_of_sequence = len(sequence)
        width_of_box = (self.available_width - 2*margin)/length_of_sequence

        x = margin
        y = y_upper - 0.15*cm

        self.canvas.setStrokeColorRGB(0,0,0, alpha = 1)
        self.canvas.setLineWidth(1)
        self.canvas.setLineCap(0)
        self.canvas.line(margin, y, self.available_width-margin, y)

        self.canvas.setLineCap(0)
        self.canvas.setLineWidth(0.2*cm)
        self.canvas.setStrokeColor(black)
        if self.CDS != '':
                CDS_str = self.CDS.split(';')

                for i in CDS_str:
                    i = i.split('-')
                    start = int(i[0])
                    end = int(i[1])
                    x_of_start = margin + width_of_box*start
                    x_of_end = margin + width_of_box*end
                    if self.reverse == 1:
                        s = x_of_start
                        x_of_start = self.available_width - x_of_end

                        x_of_end = self.available_width - s

                    self.canvas.line(x_of_start, y, x_of_end, y)


        x = margin

        for i in range(len(CDS_str)-1):
            r = int((CDS_str[i].split('-'))[1])
            l = int((CDS_str[i+1]).split('-')[0])
            if r == l:

                 x  = margin + r*width_of_box

                 if x < (self.available_width-margin) and x > margin:
                     self.canvas.setStrokeColorRGB(1,1,1, alpha=1)
                     self.canvas.setLineCap(1)
                     self.canvas.setDash(1)
                     self.canvas.setFillColorRGB(1,1,1, alpha=1)
                     self.canvas.setLineWidth(1.5)
                     self.canvas.line(x,y,x-0.06*cm, y+0.1*cm)
                     self.canvas.line(x,y,x+0.06*cm,y+0.1*cm)
                     self.canvas.line(x,y,x-0.06*cm,y-0.1*cm)
                     self.canvas.line(x,y,x+0.06*cm,y-0.1*cm)
                     self.canvas.setLineWidth(1.9)
                     self.canvas.line(x,y+0.07*cm,x, y+0.1*cm)
                     self.canvas.line(x,y-0.07*cm,x, y-0.1*cm)



'''
