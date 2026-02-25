#!/usr/bin/env python3
import argparse
import re
import cairo


# IUPAC degenerate bases from wikipedia for DNA/RNA; https://en.wikipedia.org/wiki/Nucleic_acid_notation
IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "[TU]",
    "U": "[TU]",
    "R": "[AG]",
    "Y": "[CTU]",
    "S": "[GC]",
    "W": "[ATU]",
    "K": "[GTU]",
    "M": "[AC]",
    "B": "[CGTU]",
    "D": "[AGTU]",
    "H": "[ACTU]",
    "V": "[ACG]",
    "N": "[ACGTU]",
}

# 5 colors for motifs
MOTIF_COLORS = [
    (0.70, 0.55, 0.90),  # light purple
    (0.90, 0.10, 0.10),  # red
    (0.98, 0.62, 0.12),  # orange
    (0.12, 0.56, 0.95),  # blue
    (0.15, 0.75, 0.28),  # green
]


# arguments that need to be input
def get_args():
    p = argparse.ArgumentParser(description="Draw motifs on FASTA sequences and output a png")
    p.add_argument("-f", required=True, help="FASTA file where exons are UPPERCASE and introns are lowercase")
    p.add_argument("-m", required=True, help="Motifs file with one motif per line and supports IUPAC")
    return p.parse_args()


class MotifPattern:
    """What motif we are looking for and stores regex and colors"""

    # motif_text is the original motif string (ex. "YGCY") and color is an r,g,b tuple
    def __init__(self, motif_text, color):
        self.text = motif_text.strip()
        self.color = color
        self.length = len(self.text)

        # precompile regex for this motif that is not case sensitive and allows for overlaps
        self.regex_body = self._build_regex_body()
        # lookahead so we find overlapping matches
        self.regex = re.compile(r"(?=(" + self.regex_body + r"))", re.IGNORECASE) 

    # building the regex body by converting IUPAC characters to regex and any non IUPAC characters are escaped
    def _build_regex_body(self):
        out = []
        for ch in self.text.upper():
            if ch in IUPAC:
                out.append(IUPAC[ch]) # if this character is in the IUPAC mapping, we append the corresponding regex fragment to the output list
            else:
                out.append(re.escape(ch)) # if this character is not in the IUPAC mapping, we escape it and append to the output list to support any non-IUPAC characters in motifs
        return "".join(out) # return the final regex body as a string

    # finds all start positions of this motif in a sequence where case does not matter and overlaps are allowed
    def find_hits(self, seq):
        """Return list of start positions where overlaps are allowed"""
        starts = []
        seq = seq.upper()
        matches = self.regex.finditer(seq) # find motifs in sequence

        for match in matches: # iterate over all matches and append the start position of each match to the starts list
            start = match.start()
            starts.append(start)
        return starts



class MotifHit:
    """Where a motif was found in a sequence and what motif it was, this stores motif text, start, end, color, lane"""

    def __init__(self, motif_text, start, end, color):
        self.motif = motif_text
        self.start = start
        self.end = end
        self.length = end - start
        self.color = color
        self.lane = 0  # going to assign later

class Exon:
    """class for one exon block in the sequence, with start, end, and length"""
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.length = end - start

class Intron:
    """class for one intron block in the sequence, with start, end, and length"""
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.length = end - start


# functions for reading files, finding exons and introns, assigning lanes, and drawing the image
def read_fasta(path):
    """Returns a list of (header, seq) records from multiline FASTA file and preserves case"""
    records = [] # holds (header, seq) tuples
    header = None
    pieces = [] # holds sequence pieces

    # read line by line so we can preserve case and ignore blank lines while handling multi line FASTA sequences
    for line in open(path).read().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):  # new header so we save the previous record and start a new one
            if header is not None:  # save previous record (if it exists)
                records.append((header, "".join(pieces))) # join pieces of the sequence together for the previous record
            header = line
            pieces = []
        else: # add line to current sequence pieces (case preserved)
            pieces.append(line)

    if header is not None: # save last record after loop ends
        records.append((header, "".join(pieces)))

    return records # list of (header, seq) tuples


# reads motifs from file and returns list of MotifPattern objects with colors assigned by order where it cycles through MOTIF_COLORS, if there are more motifs than colors
def read_motifs(path):
    """Return list of MotifPattern objects that are assigned colors by order"""
    motifs = []
    idx = 0 # index for cycling through MOTIF_COLORS

    for line in open(path).read().splitlines(): # read line by line so we can ignore blank lines and any # if it is present in an input motif file
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        color = MOTIF_COLORS[idx % len(MOTIF_COLORS)] # assign color by order
        motifs.append(MotifPattern(line, color)) # create MotifPattern object for this motif and add to list
        idx += 1 # increment index for next motif

    return motifs

# finds exons and introns in a sequence based on case and returns lists of Exon and Intron objects
def find_exons_introns(seq):
    """Return exons and introns as lists of Exon and Intron objects based on the case"""
    exons = []
    introns = []

    n = len(seq) # total length of the sequence
    i = 0 # index for iterating through the sequence
    while i < n: # loop through the sequence until we reach the end
        if seq[i].isupper(): # if we find an uppercase character, we are in an exon block
            start = i
            i += 1
            while i < n and seq[i].isupper():
                i += 1
            exons.append(Exon(start, i))
        else: # if we find a lowercase character, we are in an intron block
            start = i
            i += 1
            while i < n and seq[i].islower():
                i += 1
            introns.append(Intron(start, i))

    return exons, introns

# assigns lanes to motif hits so that they stack instead of overlapping and modifies hits in place and returns number of lanes used
def assign_lanes(hits):
    """Assigns lanes to motif hits so that they stack instead of overlapping. Also returns number of lanes used"""
    hits.sort(key=lambda h: (h.start, h.end)) # sort hits by start position (and end position as tiebreaker) so we can assign lanes in one pass
    lane_ends = []

    for h in hits:
        placed = False # flag to track if we have placed this hit in a lane yet
        for li in range(len(lane_ends)): # loop through existing lanes to see if we can place this hit in one of them, if it does not overlap with the last hit in that lane
            if h.start >= lane_ends[li]: # if this hit starts after the last hit in lane li ends, we can place it in this lane
                h.lane = li
                lane_ends[li] = h.end 
                placed = True # mark this hit as placed
                break
        if not placed: # if we could not place this hit in any existing lane, we need to create a new lane for it
            h.lane = len(lane_ends)
            lane_ends.append(h.end)

    return len(lane_ends)

# returns the prefix of a FASTA file path or the base name without theextension
def fasta_prefix(path):
    base = path.rsplit("/", 1)[-1]
    return base.rsplit(".", 1)[0] if "." in base else base

# draws the PNG image using cairo based on the records and motifs
def draw_png(fasta_path, records, motifs):
    # layout constants for the image that wont change for this project
    W = 1400
    margin_left = 250
    margin_right = 40
    margin_top = 90
    row_gap = 170
    exon_h = 44
    motif_h = 18
    lane_h = 8
    lane_gap = 6
    lane_offset = 18

    # legend
    legend_space = 120 + 24 * len(motifs)

    # scale based on longest gene
    max_len = 1
    for (_, seq) in records:
        if len(seq) > max_len:
            max_len = len(seq)

    usable = W - margin_left - margin_right # number of pixels we have to draw the longest sequence
    px_per_bp = usable / max_len # per base pair pixels to draw the sequences to scale

    def x(bp):
        return margin_left + bp * px_per_bp # convert a base pair position to an x coordinate in the image based on the scaling

    #  precompute per-record info
    rows = [] 
    max_lanes = 0

    for (header, seq) in records: # for each record, find the exons, introns, and motif hits and store them in a list along with the header and sequence for drawing this later
        exons, introns = find_exons_introns(seq)
        hits = [] # list to hold motif hits for this record
        for mp in motifs: # for each motif pattern, find the start positions of this motif in the sequence and create a MotifHit object for each one and add to the hits list for this record
            starts = mp.find_hits(seq)
            for st in starts:
                hits.append(MotifHit(mp.text.upper(), st, st + mp.length, mp.color)) 
        lanes_used = assign_lanes(hits) # assign lanes to motif hits
        if lanes_used > max_lanes:  
            max_lanes = lanes_used

        rows.append((header, seq, exons, introns, hits))

    lane_block = 0 # space needed for overlap lanes
    if max_lanes > 0:
        lane_block = max_lanes * lane_h + (max_lanes - 1) * lane_gap

    H = margin_top + len(rows) * row_gap + lane_block + legend_space 

    # cairo surface
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, W, int(H))
    ctx = cairo.Context(surface)

    # background to white to help with contrast
    ctx.set_source_rgb(1, 1, 1)
    ctx.paint()

    # drawing each record
    for i, (header, seq, exons, introns, hits) in enumerate(rows):
        y = margin_top + i * row_gap + 70
        name = header

        # label
        ctx.set_source_rgb(0, 0, 0) # black for the text
        ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        ctx.set_font_size(15) 
        ctx.move_to(20, y - 32) # position for the label to the left of the sequence
        ctx.show_text(name)

        # intron baseline
        ctx.set_line_width(4) 
        ctx.set_source_rgb(0, 0.0, 0.0) # black for intron
        ctx.move_to(x(0), y)
        ctx.line_to(x(len(seq)), y)
        ctx.stroke()

        # exons as filled black blocks
        ctx.set_source_rgb(0, 0, 0)
        for ex in exons:
            w = x(ex.end) - x(ex.start)
            if w < 1:
                w = 1
            ctx.rectangle(x(ex.start), y - exon_h / 2, w, exon_h) # drawing the exon as a rectangle on the base line with height exon_h and width based on the exon length that is scaled to the image
            ctx.fill()

        # motifs on main track (not offset)
        for h in hits:
            w = x(h.end) - x(h.start) # the width of the motif hit in pixels
            if w < 1:
                w = 1
            r, g, b = h.color
            ctx.set_source_rgb(r, g, b)
            ctx.rectangle(x(h.start), y - motif_h / 2, w, motif_h) # drawing the motif hit as a rectangle on the main track with the assigned color for this motif
            ctx.fill()

        # overlap lanes below
        strip_y0 = y + exon_h / 2 + lane_offset # the y coordinate of the top of the first lane strip below the exon
        for h in hits:
            lane_y = strip_y0 + h.lane * (lane_h + lane_gap) 
            w = x(h.end) - x(h.start) 
            if w < 1:
                w = 1
            r, g, b = h.color
            ctx.set_source_rgb(r, g, b)
            ctx.rectangle(x(h.start), lane_y, w, lane_h)
            ctx.fill()

    # legend 
    legend_y = H - legend_space + 40

    ctx.set_source_rgb(0, 0, 0)
    ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    ctx.set_font_size(16)
    ctx.move_to(20, legend_y - 20)
    ctx.show_text("Legend")

    ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(14)

    # intron key
    ctx.set_line_width(3)
    ctx.set_source_rgb(0.0, 0.0, 0.0)
    ctx.move_to(20, legend_y)
    ctx.line_to(60, legend_y)
    ctx.stroke()
    ctx.set_source_rgb(0, 0, 0)
    ctx.move_to(70, legend_y + 5)
    ctx.show_text("Intron")

    # exon key
    ctx.set_source_rgb(0, 0, 0)
    ctx.rectangle(20, legend_y + 15, 40, 16)
    ctx.fill()
    ctx.move_to(70, legend_y + 29)
    ctx.show_text("Exon")

    # motif legend where colors assigned by motif order
    yk = legend_y + 60
    for mp in motifs:
        r, g, b = mp.color
        ctx.set_source_rgb(r, g, b)
        ctx.rectangle(20, yk - 14, 18, 18)
        ctx.fill()

        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(50, yk)
        ctx.show_text(mp.text.upper())
        yk += 24

    out_png = fasta_prefix(fasta_path) + ".png"
    surface.write_to_png(out_png)
    print("Wrote:", out_png)

# main function to run the program
def main():
    args = get_args()
    records = read_fasta(args.f)
    motifs = read_motifs(args.m)
    draw_png(args.f, records, motifs)


if __name__ == "__main__":
    main()