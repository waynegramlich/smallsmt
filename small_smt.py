#!/usr/bin/env python

#<----------------------------------------- 100 characters --------------------------------------->|

# Classes are listed alphabetically.

from EZCAD3 import *

class Frame(Part):
    """ *Frame*: Represents the relevant *Frame* portion of pick and place machine.
    """

    def __init__(self, up, name):
        """ *Frame* Initialize the *Frame* object (i.e. *self*):
	"""

	# Standard initialization sequence:
	frame = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	Part.__init__(frame, up, name)

	# Create the various portions of the *frame*:
	frame.east_rail_ = Rail(frame, "East_Rail", True)
	frame.west_rail_ = Rail(frame, "West_Rail", False)

    def construct(self):
	""" *Frame*: Construct the *Frame* assembly object (i.e. *self*.)
	"""

	# Specify various frame constants here:
	frame = self
	frame.rail_dx_l        = rail_dx        = L(mm=9.20)
	frame.rail_dx_pitch_l  = rail_dx_pitch  = L(mm=216.0) + rail_dx
	frame.rail_dy_l        = rail_dy        = L(mm=325.0)
	frame.holes_pitch_dy_l = holes_pitch_dy = L(mm=20.00)
	frame.rail_top_z_l     = rail_top_z     = L(mm=-8.00)
	frame.holes_dy_l       = holes_dy      = float(int(rail_dy/holes_pitch_dy)) * holes_pitch_dy

class Rail(Part):
    """ *Rail*: Represents one of the rails.
    """

    def __init__(self, up, name, is_east):
	""" *Rail*: Initialize the *Rail* object (i.e. *self*.)
	"""

	# Standard initialization sequence:
	rail = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(is_east, bool)
	Part.__init__(rail, up, name)

	# Hang onto the *is_east*:
	rail.is_east_b = is_east

    def construct(self):
	""" *Rail*: Construct the *Rail* object (i.e. *self*.)
	"""

	# Grab some values from *frame*:
	rail          = self
	is_east       = rail.is_east_b
	frame         = rail.up
	rail_dx       = frame.rail_dx_l
	rail_dx_pitch = frame.rail_dx_pitch_l
	rail_dy       = frame.rail_dy_l
	rail_top_z    = frame.rail_top_z_l
	holes_dy      = frame.holes_dy_l

	# Load some constants into *rail*:
	rail.dx_l             = dx             = rail_dx
	rail.dy_l             = dy             = rail_dy
	rail.dz_l             = dz             = L(mm=7.00)
	rail.holes_pitch_dy_l = holes_pitch_dy = L(mm=20.00)

	# Compute some X coordinates:
	zero = L()
	x10 =  rail_dx_pitch/2 + dx/2
	x8  =  rail_dx_pitch/2
	x6  =  rail_dx_pitch/2 - dx/2
	x5  =  zero
	x4  = -rail_dx_pitch/2 + dx/2
	x2  = -rail_dx_pitch/2
	x0  = -rail_dx_pitch/2 - dx/2

	# Compute some Y coordinates:
	y10 =  dy/2
	y8  =  holes_dy/2
	y5  =  zero
	y2  = -holes_dy/2
	y0  = -dy/2

	# Compute some X coordinates:
	z10 = rail_top_z
	z5  = rail_top_z - dz/2
	z0  = rail_top_z - dz

	# Create the Rail as a block:
	material = Material("Plastic", "HDPE")
	color = Color("brown")
	if is_east:
	    x_west = x6
	    x_center = x8
	    x_east = x10
	else:
	    x_west = x0
	    x_center = x2
	    x_east = x4

	# Start with a block of *material*:
	corner1 = P(x_west,  y0,  z0)
	corner2 = P(x_east, y10, z10)
	rail.block("Block", material, color, corner1, corner2, "")
	rail.vice_mount("Top_Vice", "t", "n", "", zero, zero)	

	index = 0
	y = y2
	while y <= y8:
	    start = P(x_center, y, z10)
	    stop1  = P(x_center, y, z5)
	    stop2  = P(x_center, y, z0)
	    comment1 = "Hole {0}a".format(index)
	    comment2 = "Hole {0}2".format(index)
	    diameter1 = dx/2
	    diameter2 = dx/4
	    rail.hole(comment1, diameter1, start, stop1, "f")
	    rail.hole(comment2, diameter2, start, stop2, "t")
	    y += holes_pitch_dy

class SmallSMT(Part):
    """ *SmallSMT*: Represents the entir small SmallSMT pick-and-place machine of interest.
    """

    def __init__(self, up, name):
        """ *SmallSMT*: Initialize the *SmallSMT* object (i.e. *self*.)
	"""

	# Standard initialization sequence:
	small_smt = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	Part.__init__(small_smt, up, name)

	small_smt.frame_ = Frame(small_smt, "Frame")
	small_smt.trays_ = Trays(small_smt, "Trays")

    def construct(self):
        """ *SmallSMT*: Construct the *SmallSMT* object (i.e. *self8:)
	"""

	pass

class Tray(Part):
    """ *Tray*: Represents one cut-tape tray. """

    def __init__(self, up, name, colors, first_hole, holes_count, tape_widths):
	""" *Tray*: Initialize *Tray* assembly object.
	"""
	
	# Standard initialization sequence:
	tray = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	Part.__init__(tray, up, name)

	# Verify addition arguments:
	assert isinstance(colors, tuple) and len(colors) == 3
	assert isinstance(first_hole, int) and first_hole >= 0
	assert isinstance(holes_count, int) and holes_count >= 1
	assert isinstance(tape_widths, tuple) and len(tape_widths) >= 1

	# Define the tray components:
	tray.base_ = Tray_Base(tray, name + "_Base", Color(colors[0]))
	#tray.lip_  = Tray_Lid(tray,  name + "_Lid", Color(colors[1]))
	#tray.top_  = Tray_Top(tray,  name + "_Top", Color(colors[2]))

	# Save the additional arguments into *tray*:
	tray.first_hole_i  = first_hole
	tray.holes_count_i = holes_count
	tray.tape_widths_o = tape_widths

    def construct(self):
        """ *Tray*: Constructs the *Tray* object. """

	# Grab some values from *tray*:
	tray        = self
	first_hole  = tray.first_hole_i
	holes_count = tray.holes_count_i
	tape_widths = tray.tape_widths_o

	# Grap some useful *Part*'s:
	trays     = tray.up
	small_smt = trays.up
	frame     = small_smt.frame_
	west_rail = frame.west_rail_
	
	rail_dx        = frame.rail_dx_l
	rail_dx_pitch  = frame.rail_dx_pitch_l
	holes_pitch_dy = frame.holes_pitch_dy_l
	holes_dy       = frame.holes_dy_l

	zero = L()
	tray.dx_l       = dx       = rail_dx_pitch + 2 * rail_dx
	tray.dy_l       = dy       = (holes_count + 1) * holes_pitch_dy
	tray.north_y_l  = north_y  = -holes_dy/2 + (first_hole + holes_count + 0.5) * holes_pitch_dy
	tray.south_y_l  = south_y  = -holes_dy/2 + (first_hole - 0.5) * holes_pitch_dy
	tray.dz_l       = dz       = L(inch=0.500)
	tray.top_z_l    = top_z    = zero
	tray.bottom_z_l = bottom_z = -dz

	#print("tape_widths={0}".format(tape_widths))
	tray.channel_widths_o = channel_widths = [L(mm=float(width)) for width in tape_widths]
	assert isinstance(channel_widths, list)
	zero = L()
	total_channels_width = zero
	for channel_width in channel_widths:
	    #print("channel_width={0:m}".format(channel_width))
	    total_channels_width += channel_width
	#print("total_channels_width={0:m}".format(total_channels_width))
	left_over = dy - total_channels_width
	gap_dy = left_over / (len(channel_widths) + 1)
	#print("dy={0:m} left_over={1:m} gap_dy={2:m}".format(dy, left_over, gap_dy))

	tray.channel_y_centers_o = channel_y_centers = list()
	#print("south_y={0:m}".format(south_y))
	y = south_y + gap_dy
	for channel_width in channel_widths:
	    channel_y_center = y + channel_width/2
	    #print("y={0:m} channel_width={1:m} channel_y_center={2:m}".
	    #  format(y, channel_width, channel_y_center))
	    channel_y_centers.append(channel_y_center)
	    y += channel_width + gap_dy
	#print("north_y={0:m}".format(north_y))

class Tray_Base(Part):
    """ *Tray_Base*: Represents the Tray base that holes the cut tape. """

    def __init__(self, up, name, color):
	""" *Tray_Base*: Represents
	"""

	# Standard initialization sequence:
	tray_base = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(color, Color)
	Part.__init__(tray_base, up, name)
	tray_base.color_o = color

    def construct(self):
        """ *Tray_Base*: Construct the *Tray_Base* object.
	"""

	tray_base     = self
	color         = tray_base.color_o
	tray          = tray_base.up
	tray_dx       = tray.dx_l
	tray_dy       = tray.dy_l
	tray_dz       = tray.dz_l
	tray_north_y  = tray.north_y_l
	tray_south_y  = tray.south_y_l
	tray_top_z    = tray.top_z_l
	tray_bottom_z = tray.bottom_z_l
	tray_channel_y_centers = tray.channel_y_centers_o
	tray_channel_widths    = tray.channel_widths_o

	material = Material("Plastic", "HDPE")
	corner1 = P(-tray_dx/2, tray_south_y, tray_bottom_z)
	corner2 = P( tray_dx/2, tray_north_y, tray_top_z)

	tray_base.block("Block", material, color, corner1, corner2, "")
	extra_dx = L(inch="1/4")
	extra_dy = L(inch="1/4")
	tray_base.vice_mount("Top_Vice", "t", "n", "l", extra_dx, extra_dy)

	if isinstance(tray_channel_widths, list):
	    for index, channel_width in enumerate(tray_channel_widths):
		channel_y_center = tray_channel_y_centers[index]
		corner1 = P(-tray_dx/2, channel_y_center - channel_width/2, tray_top_z - tray_dz/2)
		corner2 = P( tray_dx/2, channel_y_center + channel_width/2, tray_top_z)
		radius = L(inch="1/16")
		comment = "Channel {0}".format(index)
		#print("corner1={0:m} corner2={1:m}".format(corner1, corner2))
		tray_base.simple_pocket(comment, corner1, corner2, radius, "")

class Trays(Part):
    """ *Trays*: Represents all of the trays.
    """

    def __init__(self, up, name):
	""" *Trays: Initialize the *Trays* assembly object. """

	# Standard initialization sequence:
	trays = self
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str) and not ' ' in name
	Part.__init__(trays, up, name)

	# Define each tray:
	colors1 = ("lime", "green", "dark_green")
	trays.tray1_ = Tray(trays, "Tray1", colors1, 0, 2, (4, 4, 8, 4) )
	colors2 = ("pink", "red", "dark_read")
	trays.tray2_ = Tray(trays, "Tray2", colors2, 3, 2, (4, 8, 12) )
	colors3 = ("light_blue", "blue", "dark_blue")
	trays.tray3_ = Tray(trays, "Tray3", colors3, 6, 2, (8, 12, 16) )

    def construct(self):
	""" *Trays*: Construct the *Trays* asembly object."""

	pass

def main():
    print("Hello")
    ezcad = EZCAD3(0)
    small_smt = SmallSMT(None, "SmallSMT")
    small_smt.process(ezcad)

if __name__ == "__main__":
    main()
