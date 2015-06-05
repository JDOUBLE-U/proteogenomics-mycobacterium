# coding=utf-8
# !/usr/bin/python
import os

from weblogolib import *

"""
Jeroen Merks
Jan Willem Wijnands

Instructions:
 - Install ghostscript: http://pages.cs.wisc.edu/~ghost/
    - On Windows: Add ghostscript to your system PATH
 - Install Numpy: http://sourceforge.net/projects/numpy/files/NumPy/1.9.2/
"""


class WebLogoGenerator:
    """
    ...

    :param output_folder: str
    """

    def __init__(self, output_folder):
        self.output_folder = output_folder
        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)

        # TODO fiddle around with weblogo_options for optimal results
        self.weblogo_options = LogoOptions()
        self.weblogo_options.resolution = 500
        self.weblogo_options.fineprint = "LUMC/HSLeiden"
        self.weblogo_options.logo_title = "MetAP motif"
        self.weblogo_options.title_fontsize = 6
        self.weblogo_options.number_fontsize = 6
        self.weblogo_options.fontsize = 6
        self.weblogo_options.stack_aspect_ratio = 10
        self.weblogo_options.yaxis_scale = 1
        self.weblogo_options.stack_width = std_sizes["medium"]

    def create_weblogo(self, target_proteins):
        aa_seqs = read_seq_data(open(target_proteins))

        raw_data = LogoData.from_seqs(aa_seqs)

        weblogo_format = LogoFormat(raw_data, self.weblogo_options)
        byte_array = png_formatter(raw_data, weblogo_format)
        motif_png = open(
            self.output_folder + target_proteins[
                                 target_proteins.rfind("\\"):target_proteins.rfind(".")] + "_weblogo.png",
            'wb+')
        motif_png.write(byte_array)
        motif_png.close()


"""All physical lengths are measured in points. (72 points per inch, 28.3 points per cm)

String attributes:
    o creator_text      -- Embedded as comment in figures.
    o logo_title
    o logo_label
    o unit_name         -- See std_units for options. (Default 'bits')
    o yaxis_label       -- Defaults to unit_name
    o xaxis_label
    o fineprint         -- Defaults to WebLogo name and version

Boolean attributes:
    o show_yaxis
    o show_xaxis
    o show_ends
    o show_fineprint
    o show_errorbars    -- Draw errorbars (default: False)
    o show_boxes        -- Draw boxes around stack characters (default: True)
    o debug             -- Draw extra graphics debugging information.
    o rotate_numbers    -- Draw xaxis numbers with vertical orientation?
    o scale_width       -- boolean, scale width of characters proportional to ungaps
    o pad_right         -- Make a single line logo the same width as multiline logos (default: False)

Other attributes:
    o stacks_per_line

    o yaxis_tic_interval
    o yaxis_minor_tic_ratio
    o yaxis_scale
    o xaxis_tic_interval
    o number_interval

    o shrink_fraction       -- Proportional shrinkage of characters if show_boxes is true.

    o errorbar_fraction
    o errorbar_width_fraction
    o errorbar_gray

    o resolution             -- Dots per inch (default: 96). Used for bitmapped output formats

    o default_color
    o color_scheme

    o stack_width           --
    o stack_aspect_ratio    -- Ratio of stack height to width (default: 5)

    o logo_margin           -- Default: 2 pts
    o stroke_width          -- Default: 0.5 pts
    o tic_length            -- Default: 5 pts
    o stack_margin          -- Default: 0.5 pts

    o small_fontsize        -- Small text font size in points
    o fontsize              -- Regular text font size in points
    o title_fontsize        -- Title text font size in points
    o number_fontsize       -- Font size for axis-numbers, in points.

    o text_font
    o logo_font
    o title_font

    o first_index
    o logo_start
    o logo_end

"""