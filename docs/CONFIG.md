# Configuration file parameters


Svist4get configuration file allows detailed customization of the output images. Below, comments (*placed after ; in each line*) are allowed in configuration files. Here they are used to provide short parameter descriptions.


### General options

__;[GENERAL OPTIONS]__  
__outpit_filename = svist4get_output__⠀⠀*;default filename, overridden by '-o' command line key.*  
__png_dpi = 300__⠀⠀*;DPI for bitmap png export.*  
__revcomp_transform = 0__⠀⠀*;0 or 1, similar to '-rc' command line key.*  
__hide_introns = 0__⠀⠀*;0 or 1, similar to '-hi' command line key.*  
__show_title = 1__⠀⠀*;1 or 0, show or hide image title.*  
__font_size_title = 7.5__⠀⠀ *;image title text font size.*  
__show_aa_seq_track = auto__⠀⠀*;can be 'auto' (show or hide track based on a number of transcripts), 1 or 0, show or hide aminoacid sequence track. The track is mostly useful for studying ORFs in a particular transcript.*   
__show_genomic_axis = 1__⠀⠀*;1 or 0, show or hide the genomic coordinates axis.*  
__show_genomiс_axis_tics = 0__⠀⠀*;0 or 1, show or hide the genomic coordinates axis tics.*  
__show_nt_seq_track = auto__⠀⠀*;can be 'auto' (show or hide track based on genomic window width), 1, or 0 (force show or hide).*  
__transcript_label_style = auto__⠀⠀*;can be 'none' (hide the transcript label), 'name', 'id', 'gene_id' (show annotated transcript_name, transcript_id or gene_id respectively as a transcript_label), 'both' (show annotated transcript_name and transcript_id) and 'both' (hides name if it is included the ID)*


### File paths

**NOTE:** The file paths are given relative to the CONFIG file path.

__;[GENERAL FILE PATHS]__  
 __triplet_code = ./triplet_code.txt__⠀⠀*;the triplet code table file.*  
__palette = ./palette.txt__⠀⠀*;the color palette file, colors are defined by standard hex codes.*  
__gtf_file = -__⠀⠀*;default gff file path, overridden by '-gtf' command line key.*  
__fasta_file = -__⠀⠀*;default genome fasta file path, overridden by '-fa' command line key.*

__;[FONT FILE PATHS]__  
__mono_font = fonts/iosevka-regular.ttf__⠀⠀*;default monotype font for nucleotide and aminoacid labels.*  
__regular_font = fonts/Lato-Regular.ttf__⠀⠀*;default font for other labels.*


### Image size and style customization

**NOTE:** The page size parameters values are in centimeters. Page width determines the image width taking into account the defined margins. Page height is unlimited (i.e. it is automatically extended to fit all user-defined tracks).

__;[GENERAL OUTPUT PAGE PARAMETERS [cm]]__  
__page_width = 8__  
__margin = 0.5__  
__gap = 0.05__⠀⠀*;vertical spacer between tracks.*


__;[VERTICAL GRID SETTINGS]__  
__vgrid_step = 0.5__⠀⠀*;distance between vertical grid lines [cm].*  
__c_vgrid = blue__⠀⠀*;vertical grid color.*  
__c_vgrid_alpha = 0__⠀⠀*;vertical grid opacity, from 0 to 1. By default in the latest version it is set as 0 to hide the vertical grid.*

__;[FRAME HIGHLIGHT SETTINGS]__  
__c_fill_hframe = blue__⠀⠀*;fill color.*  
__c_fill_hframe_alpha = 0.25__⠀⠀*;fill opacity, from 0 to 1.*  
__c_stroke_hframe = blue__⠀⠀*;border color.*  
__c_stroke_hframe_alpha = 0.05__⠀⠀*;border opacity, from 0 to 1.*  
__c_text_hframe = black__⠀⠀*;text label color.*  
__hframe_font_size = 5__⠀⠀*;text label font size.*

__;[GENOMIC COORDINATES AXIS SETTINGS]__  
__gaxis_label_font_size = 7__  
__gaxis_tics_font_size = 6__  
__gaxis_tics_step = 0__⠀⠀*;the value must be in nucleotides, 0 means automatic selection.*  

__;[TRANSCRIPT STRUCTURE TRACK SETTING]__  
__stroke_width_CDS = 0.14__⠀⠀*;[cm]*  
__stroke_width_exon = 0.08__⠀⠀*;[cm]*  
__transcript_id_font_size = 5.5__  
__arrow_height = 0.03__⠀⠀*;height of arrows showing transcript orientation [cm].*  
__c_marks = deepred__⠀⠀*;color of marks used to pinpoint exon junctions when intronic regions are hidden.*


__;[BEDGRAPH TRACK SETTINGS]__  
__bedgraph_track_height = 0.75__⠀⠀*;[cm], paired bedGraph tracks will use  the same height.*  
__bedgraph_column_min_width = 0.04__⠀⠀*;[cm], minimal visible width of a single bar of the bedgraph barchart, affects automatic selection of the genomic window size that is used to aggregate values for each bar of the bedgraph track.*  
__bedgraph_tics_font_size = 5__  
__bedgraph_label_font_size = 7__  
__bedgraph_label_position = center__ *; can be 'right', 'left', or 'center'. Overridden by '-blp' command line key. *  
__c_bedgraph_label_alpha = 0.65__⠀⠀*;label opacity, from 0 to 1
c_bedgraph_tracks = pink,blue,green,yellow,purple,orange,brown,red; cycling default palette for bedgraph tracks.*  
__c_bedgraph_alpha = 0.65__⠀⠀*;bedgraph track opacity, from 0 to 1.*  
__bedgraph_bar = mean__⠀⠀*;By default, the bedGraph tracks are plotted with 'bars' spanning small genomic windows of automatically optimized length. This parameter sets the function to aggregate bedGraph signal values for a single bar of the bedGraph tracks.  Default: 'mean'. Alternatives: 'max', 'min', 'median' and 'none'. 'none' forces plotting the raw signal at 1nt resolution.*  
__bedgraph_axis_tics = auto__⠀⠀*;can be either 'auto' (drawing tics at the top and the bottom) or a complete comma-separated list of values (in the bedGraph signal values scale). Overridden by '-lb' command line parameter. Not available for paired bedGraph tracks.*  
__bedgraph_axis_tics_step = 0__⠀⠀*;can be set either to 0 (auto) or to the actual value using the same scale as that of the signal track (i.e. the bedGraph file contents)*  
__bedgraph_upper_limit = auto__⠀⠀*;can be 'auto','min', or 'max' or a particular numeric value. Overridden by '-bul' command line key. 'min' and 'max' will set a single limit for all tracks based on minimal or maximal in-window value of all bedGraph tracks. For a paired bedGraph track, this sets limit for the first track in a pair.*  
__bedgraph_lower_limit = auto__⠀⠀*;can be 'auto','min', or 'max' or a particular numeric value. Overridden by '-bll' command line key. 'min' and 'max' will set a single limit for all tracks based on minimal or maximal in-window value of all bedGraph signals. For a paired bedGraph track, this sets limit for the second track in a pair.*


__;[NUCLEOTIDE SEQUENCE TRACK SETTINGS]__  
__c_A = blue__  
__c_A_alpha = 0.8__  
__c_G = orange__  
__c_G_alpha = 0.8__  
__c_C = yellow__  
__c_C_alpha = 0.8__  
__c_T = purple__  
__c_T_alpha = 0.8__

__;[AMINOACID SEQUENCE TRACK SETTINGS]__  
__aa_line_height = 0.14__⠀⠀*;[cm], height of start and stop codon tics for the track in 'tics' style.*  
__aa_line_width = 0.05__⠀⠀*;[cm], width of tics for 'tics' style.*  
__aa_show_start = 1__⠀⠀*;1 or 0, show or hide start codons.*  
__aa_track_style = auto__⠀⠀*;can be 'auto', 'codons', or 'tics'. 'auto' selects the representation based on the visible width of the genomic region.*  
__c_stop_codon = deepred__  
__c_stop_codon_alpha = 0.59__  
__c_start_codon = deepgreen__  
__c_start_codon_alpha = 0.59__  
__c_highlight_codon = blue__⠀⠀*;color for user-selected codons 'highlighted' with '-hc' command line key.*  
__c_highlight_codon_alpha = 0.7__⠀⠀*;opacity for user-selected codons 'highlighted' with '-hc' command line key.*  
__c_other_codon = lightgray__  
__c_other_codon_alpha = 0.7__  
__c_part_codon = gray__  
__c_part_codon_alpha = 0.8__  
__c_aa_track_tics_line = lightgray__⠀⠀*;color of the track line for the 'tics' style.*  
__c_aa_track_tics_line_alpha = 0.6__  
__c_fill_aa_frame = blue__  
__c_fill_aa_frame_alpha = 0.2__  
__c_stroke_aa_frame = blue__  
__c_stroke_aa_frame_alpha = 0.6__

__;[GENOMIC REGIONS TRACK SETTINGS]__  
__regions_line_width = 0.1__  
__font_size_regions_label = 6__  
__c_regions = black__⠀⠀*;it is possible to provide a comma-separated list of colors in a way similar to that of bedgraph tracks.*  
__c_regions_alpha = 0.55__



### [Main page](https://github.com/art-egorov/svist4get)  
#### [Quickstart guide](./QSGUIDE.md)  
#### [Command-line parameters](./PARAMETERS.md)  
#### [API usage examples](./API.md)  
#### [FAQ](./FAQ.md)  
#### [Version log](./VERSION.md)

